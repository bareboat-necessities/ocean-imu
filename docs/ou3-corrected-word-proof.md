# Corrected-word covariance reduction and finite-error composition

Regime scope: `ou3-regime-design.md` separates stationary observability and
finite transitions from complete excited moving windows on this same carried
proof path. `ou3-moving-six-pivots.md` supplies a conditional all-six-pivot
bound with chronological reset and actual-row defects. Its uniform premises
remain OPEN; neither this reduction nor the mode split certifies B_*, J_AG,
rho_0, a nonlinear retained radius or every-prefix retention.


These results enter the single tail inequality through

`sqrt(V_end) <= sqrt(1-delta) sqrt(V_root) + sum_i s_i`.

They retain the complete 21-state covariance, actual gains, interleaved
corrections and resets. They do not yet supply a uniform `delta`, capture or
an invariant nonlinear region. All statements below concern regular default
real-arithmetic A21, with no frame/relock reconfiguration. Float32 is separate.

## 1. An embedded nuisance floor at the correct word root

Partition `x=(h,n)`, where `h=(theta,b_g)` has six coordinates and
`n=(v,p,S,a_w,b_a)` has fifteen. Start each word immediately **before** a
prediction, after at least 17 s of regular A21. This is a choice of where to
evaluate the carried history, not a covariance or state restart.

The existing joint comparison, immediately after the preceding prediction,
implies `P >= E_n L_post E_n'`. In raw physical coordinates its LIN part is
`D A^-1 D/2`, with `D=diag(2.4,18,132,4)` and the retained full 4x4 action
matrix A, tensored with I3. Its BA part is `q_BA I3/2`.

For a correction, the covariance map

`C_H(P)=P-P H'(H P H'+R)^-1 H P`

is monotone on positive semidefinite matrices. This follows by minimizing
the covariance of the linear estimation error over the gain; equivalently,
apply inverse order for positive definite P and take a positive-definite
limit. In particular, evaluating at a singular embedded nuisance covariance
gives

`C_H(P) >= E_n (L^-1+H_n' R^-1 H_n)^-1 E_n'`.

The actual accelerometer nuisance row is `[0,0,0,R_wb,I]`, so

`H_n' R_acc^-1 H_n <= (2/.05^2) diag(0,0,0,I3,I3)`.

The S row adds at most `I3/.075^2` on S. There is at most one of each between
successive predictions. Magnetic rows have `H_n=0`. Every literal attitude
reset has `G E_n=E_n`; arbitrary such resets therefore preserve this embedded
floor. Actual PSD synchronization increments can only increase it. No bound
on the nominal AW mean, attitude injection or AG covariance is used here.

Consequently the nuisance floor at the next pre-prediction root is

`L = (L_post^-1 + diag(0,0,I3/.075^2,2 I3/.05^2,2 I3/.05^2))^-1`.

The existing nuisance upper comparison supplies `P_nn <= U` at this same
root. Both L and U are positive definite, uniform matrices. The exact
certificate supplies `0<alpha<1` with `L >= alpha U`. It preserves the full
4x4 LIN precision and verifies the matrix inequality rationally; it does not
interpret LDL pivots as eigenvalues.

## 2. Six corrected loss columns suffice with these covariance bounds

**Lemma.** Suppose

`P > 0`, `P >= E_n L E_n'`, `P_nn <= U`, `L >= alpha U`, `0<alpha<1`,

and the **complete corrected word loss** satisfies

`D_word,hh >= J > 0`.

Then, for every eta>0,

`P <= C_eta := diag((1+eta) J^-1/alpha, (1+1/eta) U)`.

**Proof.** Write `P=[[A,C],[C',N]]`. Since `L >= alpha U >= alpha N`,

`[[A,C],[C',(1-alpha)N]] >= 0`.

Taking its Schur complement yields

`C N^-1 C' <= (1-alpha) A`,

and hence the actual conditional covariance obeys

`S_h := A-C N^-1 C' >= alpha A`.

The complete energy identity gives `0 <= D_word <= P^-1`. Its principal
six-coordinate restriction therefore gives

`S_h^-1 = (P^-1)_hh >= D_word,hh >= J`.

Thus `A <= J^-1/alpha`. For an arbitrary full vector, the PSD cross-block
Cauchy inequality and Young inequality give

`[h;n]'P[h;n] <= (1+eta)h'Ah + (1+1/eta)n'Nn`.

Substitution proves the claimed full upper comparison. No principal loss is
embedded as an independent measurement, and C is never set to zero.

At the first prediction of this word, suppose its actual process covariance
satisfies

`Q >= epsilon F C_eta F'`, `epsilon>0`.

Then `P_minus >= (1+epsilon) F P F'`, so the first prediction alone gives

`D_word >= [epsilon/(1+epsilon)] P^-1`.

Every subsequent regular correction/reset/prediction preserves the
nonincrease of the homogeneous covariance-weighted energy. Therefore

`rho_0 <= 1/(1+epsilon) < 1`.

This is a full 21-state implication. It applies uniformly **if J is proved
uniformly** for the six AG columns of the actual transported loss. The
two-column MAGNETIC SERVICE restriction does not supply J. Neither does the
physical three-dimensional vector Gramian. Those are still distinct objects.

For existence, the source has a uniform positive full Q floor and a uniform
bound on F. In the LIN block, with `x=h/tau` and natural step scales
`D_h=diag(h,h^2,h^3,1)`,

`Q_ideal = sigma^2 x D_h B(x) D_h`, `B(x) >= B(0)/(1+x)^2`.

The latter comparison follows by expressing the undamped impulse polynomial
as `(I+x V)` times the damped impulse, where the Volterra integration
operator on [0,1] has norm at most one. Charge the already bounded source
small-x covariance defect. Fresh AG and BA process floors complete all 21
coordinates. The literal transition obeys `||F|| <= 1.012` for
`h in [.004,.006]`; the exact attitude rotation has norm one and its gyro-bias integral has norm
at most h. In the small-rate polynomial branch, bound these by 1+h/2 and
3h/2, respectively; their sum is at most 1+2h. The LIN row and
column sums are at most `1+h+h^2/2+h^3/6`, and BA decay is at most one.

The scalar Q floor is an existence bound, not a practically useful nonlinear
margin. For an explicit margin retain the block matrix inequality for epsilon.
The accompanying 80-digit feasibility experiment and independent exact
cross-coupled example pass; neither is a shipping-uniform J certificate.

## 3. Actual-gain finite-error composition

On one realized estimator execution, freeze its actual coefficients only for
the auxiliary linear error comparison. The true finite error satisfies

`e_i = A_i e_(i-1) + d_i`.

No second estimator is run, and no equality of two estimators' gains or
event decisions is assumed. A_i includes the realized prediction, optimal
Joseph correction with its effective R, or literal reset. Each obeys

`A_i' P_i^-1 A_i <= P_(i-1)^-1`.

Variation of constants and the triangle inequality consequently give

`||e_end||_(P_end^-1) <= ||M e_root||_(P_end^-1)
                       + sum_i ||d_i||_(P_i^-1)`.

If the complete loss has margin delta, the first term is at most
`sqrt(1-delta)||e_root||_(P_root^-1)`.
The same proof at every prefix uses coefficient one in place of the strict
word coefficient. This is an exact finite-error statement whenever each d_i
is the difference from the actual operation, including projection.

For a correction with innovation residual `H e + r`, the additive defect
before reset is `-K r`. The Joseph covariance satisfies `P_plus >= K R K'`,
therefore

`||K r||_(P_plus^-1) <= ||r||_(R^-1)`.

This avoids a separate gain-norm bound or a derivative of the gain. R must
include the actual innovation safety increment. For the lever-arm-disabled
accelerometer, `f_hat=R_hat(a_hat-g)` and the left attitude error theta give

`||r_acc,nonlinear|| <= ||f_hat|| |theta|^2/2 + |theta| |e_aw|`.

This follows from `||exp([theta])-I-[theta]|| <= |theta|^2/2` and
`||exp([theta])-I|| <= |theta|`. The analogous magnetic curvature term is
`|B_hat| |theta|^2/2`; reference error and physical magnetic residuals remain
additional forcing. Prediction mismatch uses the established joint metric
action bound. Euclidean BA and BG projections retain separate sectors/defects;
it is not declared nonexpansive in the full covariance metric.

There is a stronger accumulation bound for prediction and correction inputs.
Write the exact comparison covariance recursion as

`P_i=A_i P_(i-1) A_i' + B_i B_i'`,

where B_i is a factor of Q_i at prediction, `K_i R_i^(1/2)` at correction,
and zero at a congruent reset. Let the corresponding deterministic finite
defect be `B_i u_i`; thus its action is `d_i' Q_i^-1 d_i` at prediction and
`r_i' R_i^-1 r_i` at correction (use the least-norm factor input if singular).
Expanding the recursion through the whole word gives the matrix identity

`P_W = M P_0 M' + sum_i M_(W<-i) B_i B_i' M_(W<-i)'`.

After whitening by the terminal covariance factor, the horizontally stacked
input matrix has norm at most one. Hence the sharper finite-error inequality
is

`sqrt(V_W) <= sqrt(1-delta) sqrt(V_0) + sqrt(sum_i |u_i|^2)
             + sum_(reset/projection i) ||d_i||_(P_i^-1)`.

This is deterministic matrix algebra: the physical input terms need not be
independent or obey Gaussian/OU laws. All correlations in Q, R and P remain.
It replaces an unnecessary factor proportional to the square root of the
number of operations in the bound obtained by summing individual prediction
and correction norms. The actual inputs, which depend on the realized finite
error, must still satisfy the action bounds on the retained region. The
finite-angle reset and mean-only projection defects do not have Joseph noise
channels; they are not included in this square-summed term.

## 4. Finite-angle remainder of the literal reset

Let d be the actual requested attitude injection and `v=theta-d`, with
`|v|<=r<2`. For ideal exponential injection define

`u=Log(exp([theta]) exp(-[d]))`, `G(d)=I+[d]/2`.

Then

`|u-G(d)v| <= |d|^2 |v|/6 + (1/[2(2-r)]+1/4)|v|^2`.

To prove this, follow
`u(t)=Log(exp([d+t v]) exp(-[d]))`. Its rotation angle is at most t|v|,
since the left Jacobian `J_l(z)=int_0^1 exp(s[z]) ds` has norm at most one.
The path remains in the principal log chart and

`u'(t)=J_l(u(t))^-1 J_l(d+t v) v`.

The integral representation gives
`||J_l(z)-J_l(w)|| <= |z-w|/2` and `||J_l(u)-I|| <= |u|/2`.
The inverse Neumann bound yields

`||J_l(u(t))^-1 J_l(d+t v)-J_l(d)||
 <= t|v| (1/(2-r)+1/2)`.

Integrate and add `||J_l(d)-G(d)|| <= |d|^2/6`.
In particular, the reset remainder contains a term proportional to
`|d|^2 |v|`: with a nonzero realized injection it cannot simply be called
quadratic in the estimation error with no injection dependence.

For `|d|<.01`, the shipping quaternion uses normalized truncated polynomials,
not the exact exponential. Its rotation-angle defect is bounded by

`epsilon_poly <= 4 (|d|^6/46080 + |d|^7/645120)`.

The alternating sine/cosine remainders bound the raw quaternion error by the
parenthesized quantity. The segment from that quaternion to the unit exact
quaternion stays within 1/2 of the unit sphere, so the derivative of its
rotation angle is bounded by four. Normalizing does not change that angle.
If `r+epsilon_poly<2`, add

`2 epsilon_poly/(2-r-epsilon_poly)`

to the reset remainder above. This follows by transporting the injection
defect along its shortest rotation path and bounding the log derivative.
The real sine/cosine branch has no such polynomial defect. Floating-point
evaluation, normalization and covariance arithmetic remain separate supplies.

## 5. What remains to reach capture and retention

The reduction closes the **implication**, not its six-column J premise.
Next bound those columns under the admitted varying nominal dynamics and
actual resets. Keep their nuisance coupling through the proven covariance
comparison; do not reintroduce independent heading information.

If the summed word supplies are bounded by `a sqrt(V_root)+b`, define
`q=sqrt(1-delta)+a`. An invariant root radius requires `q<1` and
`b <= (1-q)r_root`. Every-prefix supplies must also keep all operations in
the region in which their bounds were established. Actual startup, finite
H18 and release must enter that region. None of these numerical inequalities
has yet been discharged for the declared disturbance envelopes.

The added gyro mean projection has radius .5 rad/s and physical radius .02,
so its ideal Euclidean sector gap is .48 rad/s. Carry its signed defect and
inward-rounding supply in the same sum. No inherited AG covariance ceiling or
all-prefix gyro-projection inactivity is inferred from the BA guard below.
`ou-gyro-bias-projection.md` supplies the one-prediction transport result and
states why the full signed temporal margin is still open.

## 6. Word contraction as one information-ratio inequality

Role in the tail inequality: the lemmas below bound the linear factor `rho_0`
of a complete A21 word in `V_(j+1) <= rho V_j + c_d ||d||^2`. They keep the
root covariance, cross covariance, actual gains and resets of that word.

**Smoother identity (exact).** Freeze the word's coefficients and use the
optimal gains. Let `C_k = Cov(x_0, x_k | y)` be the fixed-point smoother
cross covariance. It evolves as `C<-C F'`, `C<-C (I-KH)'` and `C<-C G'`, and
`Sigma_00` loses `C H' S^-1 H C'` at each correction. The closed-loop word
map `M` then satisfies

`M = C_end' P_0^-1`,
`M' P_end^-1 M = P_0^-1 (Sigma_00|y - Sigma_00|y,x_end) P_0^-1`.

Hence `rho_W = lambda_max(P_0^-1/2 (Sigma_00|y - Sigma_00|y,x_end) P_0^-1/2)`.
Write `J` for the root information carried by the word's data and `A` for
the root information in the data together with the terminal state. Then
`Sigma_00|y = (P_0^-1+J)^-1`, `Sigma_00|y,x_end = (P_0^-1+A)^-1` and `J<=A`.
The unknowns `w, v` are independent of `x_0`; realized coefficients are
frozen and no stochastic law is imposed on physical histories.
`word_smoother_identity` checks both identities exactly on the supplied
21-state word.

Two separated sufficient conditions follow at once:

- **Information:** `Sigma_00|y <= (1-delta) P_0`, i.e. `J >= delta/(1-delta) P_0^-1`.
- **Forgetting:** `Sigma_00|y,x_end >= delta P_0`, i.e. `A` small against `P_0^-1`.
  The first-prediction route of §2 is the special case with `J` dropped and
  one-step hiding `A <= (F^-1 Q F^-T)^-1`.

**Neither separated condition suffices (exact and carried).**
`mixed_mechanism_example` builds a two-state word with `rho=1/2` in which
both separated margins are zero. The carried shipping words show the same
split, in the replay setting of `ou3-ag-readout-proof.md`:

| Word | Quiet exact | Quiet info-only | Quiet forget-only | Wave exact | Wave info-only | Wave forget-only |
|---:|---:|---:|---:|---:|---:|---:|
| 0.32 s | 4.06e-4 | 3e-26 | 2.0e-10 | 6.98e-4 | 3.0e-10 | 2.1e-9 |
| 4 s | 0.0242 | 3e-17 | 6.3e-4 | 0.0248 | 5.5e-7 | 2.9e-3 |
| 16 s | 0.0982 | 1e-14 | 5.3e-3 | 0.0996 | 3.2e-4 | 6.6e-3 |
| 64 s | 0.389 | 3e-12 | 5.2e-3 | 0.393 | 3.0e-3 | 6.5e-3 |

The slowest 0.32-s direction is translational: about 87% velocity/position
with 10–13% accelerometer bias. Its loss comes mostly from the S
pseudo-observation and the AW sync. Bias directions contract about 3–5×
faster. Quiet water keeps the tilt/BA kernel about the magnetic axis
unobservable, so its information-only margin stays at the rounding floor.

**Information-ratio lemma (exact).** Let `C >= P_0` be any root covariance
upper bound. For `kappa >= 1` put

`k(kappa) = max(0, lambda_max(C^1/2 (A - kappa J) C^1/2))`.

Then

`rho_W <= sup_(y>=0) [ 1/(1+y) - 1/((1+k)(1+kappa y)) ] < 1`.

If `A <= kappa J` (then `k=0`), this gives `rho_W <= (sqrt(kappa)-1)/(sqrt(kappa)+1)`
for every `P_0`, with no covariance bound at all.

*Proof.* Put `Y=P_0^1/2 J P_0^1/2` and `Z=P_0^1/2 A P_0^1/2`, so that
`rho_W=lambda_max((I+Y)^-1-(I+Z)^-1)`. Let `K` be the positive part of
`A-kappa J` in the `C` metric, so `A <= kappa J + K`. Since `P_0 <= C`,
`P_0^1/2 K P_0^1/2 <= k I`. Therefore
`I+Z <= (1+k)(I+kappa Y)`, and operator monotonicity of inversion gives
`(I+Z)^-1 >= (I+kappa Y)^-1/(1+k)`. The remaining matrix is a function of
`Y` alone, so functional calculus gives the scalar supremum. With `k=0` its
maximum is at `y=1/sqrt(kappa)`. ∎

**How the pieces fit.** The lemma uses both mechanisms at once. Directions
the data see contract through `kappa`, with no covariance bound. Directions
the data miss contract through the upper comparison `C`, which enters only
in `k`. The joint reader supplies that `C` (`ou3-ag-readout-proof.md`). On
carried words its full bound is within 55–77× of the actual root covariance
once the preceding window is 16 s or longer. The first-prediction chain is
the special case `J=0`, `A <= (F^-1QF^-T)^-1`, giving `1/(1+epsilon)`.

The source-uniform obligations of this lemma are restated in section 7,
where they reduce to one diameter of the word itself.

## 7. Word Riccati diameter and the kernel-bounded contraction

Role in the tail inequality: these results bound the linear factor `rho_0`
of one complete word in `V_(j+1) <= rho V_j + c_d ||d||^2`. They retain the
actual frozen coefficients, gains, resets and correlated sources, and need at
most one scalar covariance bound. Exact checks are in `word_diameter.py`
(`word-diameter-certificate.json`); the carried feasibility test is
`information_ratio_source_diagnostic.py`
(`information-ratio-source-feasibility.json`).

**Theorem D (diameter, exact).** Write `Pi = Ric_W(0)` for the terminal
covariance of the word started from a known root and
`P_diff = lim_t Ric_W(t I)` for the diffusely started word. With
`Phi_tilde = Phi - T A' Sigma^-1 O` (the root coefficient of `E[x_N|y,x_0]`):

- `A = J + Phi_tilde' Pi^-1 Phi_tilde`, `Pi = T (I - A' Sigma^-1 A) T'`;
- `Ric_W(t I) = Pi + Phi_tilde (I/t + J)^-1 Phi_tilde'`, so
  `P_diff = Pi + Phi_tilde J^-1 Phi_tilde'` when `J > 0`;
- `kappa_W := lambda_max(J^-1 A) = lambda_max(Pi^-1 P_diff)`;
- for every root covariance, `Pi <= P_end <= P_diff` and
  `rho_W <= (sqrt(kappa_W)-1)/(sqrt(kappa_W)+1) = tanh(log(kappa_W)/4)`.

*Proof.* The first two items are the Schur and Woodbury forms of the joint
Gaussian `(y, x_N)`. With `X = Pi^-1/2 Phi_tilde J^-1/2`, the matrices
`X'X = J^-1/2 (A-J) J^-1/2` and `X X' = Pi^-1/2 (P_diff-Pi) Pi^-1/2` share
their nonzero spectrum, so both largest eigenvalues equal `kappa_W - 1`.
The bound is the section 6 lemma with `k = 0`. ∎

No covariance ceiling enters. The bound has the Birkhoff–Hopf form
`tanh(Delta/4)` at `Delta = log kappa_W`, and it is sharp. For the scalar word
"correct `y=x+v`, `R=1`, then predict with `Q=1`", `kappa=2` and
`rho(p)=p/((p+1)(2p+1))`. The identity
`(3-2 sqrt2)(2p^2+3p+1) - p = ((2-sqrt2)p - (sqrt2-1))^2` shows
`sup_p rho = (sqrt2-1)/(sqrt2+1)`, attained at `p = 1/sqrt2`.

The section 6 supremum has a closed form. With `c = 1/(1+k)` it equals
`(sqrt(kappa)-sqrt(c))^2/(kappa-1)` if `c kappa >= 1`, and `1-c` otherwise.

**Composition and factorization (exact).**

- *Composition.* For `W = W_2 o W_1`, `[Pi_W, P_diff,W]` lies inside
  `[Pi_W2, P_diff,W2]` by Riccati monotonicity. Hence `kappa_W <= kappa_W2`,
  and along an execution `rho_W <= rho_W1 rho_W2`.
- *Factorization.* Let `A_s` add the slow terminal state `(theta, b_g, b_a)`
  to the data, so `J <= A_s <= A`. Then `kappa_W <= kappa_s kappa_f|s`.
  Here `kappa_s` is the diameter of the slow terminal marginals of
  `[Pi, P_diff]` and `kappa_f|s` that of the fast Schur complements. The
  reverse order holds as well.
- *Dual certificates.* Every compression `L` of the data gives
  `J >= (L O)' (L Sigma L')^-1 (L O)`. Every full-root reader with
  `L O = Phi` gives `P_diff <= (T - L A)(T - L A)'`, with equality at the
  minimum-action reader.

**Corollary K (rank-structured kernel, exact).** Let `nu` be a root direction
and suppose `A <= kappa J + lambda nu nu'`. Then every root with
`nu' P_0 nu <= c` satisfies

`rho_W <= sup_y [ 1/(1+y) - 1/((1+lambda c)(1+kappa y)) ]`.

This is section 6 with `K = lambda nu nu'`, whose `P_0` congruence has the
single eigenvalue `lambda nu' P_0 nu`. Only the scalar `c` bounds the root
covariance; `J` may be singular along `nu`. For `mu > 0` let
`P_nu = Ric_W(root information mu nu nu') = Pi + Phi_tilde (J + mu nu nu')^-1 Phi_tilde'`
and `kappa_nu = lambda_max(Pi^-1 P_nu)`.

- *Premise supplied by a diameter.* A fictitious root observation of
  `nu'x_0` leaves `Pi` unchanged, so Theorem D gives
  `A - kappa_nu J <= (kappa_nu - 1) mu nu nu'`. With `mu = 1/c`,
  `rho_W <= 1 - 1/kappa_nu`.
- *Invariance.* If `nu' P_0 nu <= c`, then `P_0^-1 >= nu nu'/c`, so
  `P_end <= P_nu` (Riccati monotonicity, with a limit in `c`). Along a chain
  of words with kernels `nu_j` and ceilings `c_j`, the property
  `nu_j' P_j nu_j <= c_j` of the actual roots therefore propagates whenever
  `nu_(j+1)' P_nu_j nu_(j+1) <= c_(j+1)`.

`kernel_bounded_certificate` verifies the corollary exactly on a word whose
`J` is singular along `nu`.

**Physical kernel and scalar ceiling.** In A21 take
`nu = (theta_hat, 0, 0, 0, 0, 0, nu_ba)`, with `theta_hat` the unit body field
axis (the null direction of the magnetic attitude Jacobian) and
`nu_ba = -J_att theta_hat`. This is the quiet tilt/accelerometer-bias
direction that neither the field nor the accelerometer row sees. Block
Cauchy–Schwarz and the proved `P_ba,ba <= I/1600` give, for every
covariance,

`nu' P nu <= ( sqrt(tau_theta) + |nu_ba|/40 )^2`, `tau_theta >= theta_hat' P_tt theta_hat`.

With `|nu_ba|` near `g`, the BA term dominates. A tilt premise
`tau_theta = 10^-3 rad^2` raises the ceiling by 27% at `|nu_ba| = g`.

**S-chain cancellation (exact).** Use the literal per-axis LIN transition
with any coefficients `phi_va, phi_pa, phi_Sa, alpha`, fresh noise
`(n_v, n_p, n_S, n_a)` and PSD sync jumps on `a`. Let `c` annihilate
`(1, t, t^2)` on the applied S times, and put

`omega_k = sum_(t_j >= t_(k+1)) c_j (phi_Sa + phi_pa d + phi_va d^2/2)`,
`g_k = sum_(t_j >= t_(k+1)) c_j (d^2/2, d, 1)`, with `d = t_j - t_(k+1)`.

Then exactly

`sum_k omega_k a_k - sum_j c_j S_j = - sum_k g_k' (n_v, n_p, n_S)_k`.

The whole `(v, p, S, a_w)` root, every AW noise and every sync jump cancel.
With `J_aw = R_wb`, the world-frame accelerometer combination
`sum_k omega_k R_k' y_acc,k - sum_j c_j y_S,j` therefore sees only the
attitude, BA, measurement noise and this residual. Its variance is at most
`(1+eps) q sum_k |g_k|^2 (h^3/3 + h^5/20 + h^7/252)`, `q <= 2 sigma^2/tau`.

**Gyro-bias persistence cap (source-uniform).** Every regular A21 word of
duration `T` (default profile, zero lever arm) satisfies

`kappa_W >= sigma_g^2 (1-delta) / ((1+e)^2 b0 T^2)`.

Here `e <= 2/22!` bounds the degree-18 series for `B_step` and `int B`, and
`delta = h^2 b0 (1+e)^2/(4 sigma_g^2)`.

*Proof.* Take a root gyro-bias perturbation `beta`. The data-only mimic
keeps `b_g = 0` and feeds the literal attitude noise `B_step beta` at each
prediction. It reproduces every row, reset and nuisance state exactly; the
Schur complement `Q_tt - Q_tb Q_bb^-1 Q_bt >= sigma_g^2 h (1-delta)` then
gives `u'J u <= T(1+e)^2|beta|^2/(sigma_g^2(1-delta))`. The terminal bias
alone has `A >= E_bg E_bg'/(b0 T)`. ∎

The cap is 71.19, 17.80 and 4.449 at 16, 32 and 64 s. No `k = 0`
certificate on such a word can exceed the margin `2/(1+sqrt(cap))`: .212,
.383 and .643.

**Carried feasibility (non-promoting).** Float64 optimal-gain replays of
the literal exported coefficients record a 64-s history from 225 s or
5000 s, so every word starts at a root at 289 s or 5064 s
(`information-ratio-source-feasibility.json`). BA release is at 120 s. At
289 s the BA marginal (5.5e-5) is still an order of magnitude below its
5064-s level (5.3e-4 quiet, 2.5e-4 wave). Margins are `1-rho`.
"Kernel, exact" uses the actual `nu'P_0 nu` in Corollary K; "kernel,
ceiling" uses `(sqrt(10^-3)+|nu_ba|/40)^2`.

| Word | Exact | Ideal `C=P_0` | Joint-reader `C` (loss) | `kappa_W` | Kernel, exact | Kernel, ceiling |
|---|---:|---:|---:|---:|---:|---:|
| quiet 289 s, 16 s | .0676 | .0674 | .0013 (51.7×) | ∞ | .0357 | .0051 |
| quiet 289 s, 64 s | .2265 | .2262 | .0053 (43.1×) | ∞ | .2196 | .0203 |
| quiet 5064 s, 16 s | .0075 | .0075 | .0013 (5.7×) | ∞ | .0074 | .0051 |
| quiet 5064 s, 64 s | .0296 | .0295 | .0053 (5.6×) | ∞ | .0295 | .0203 |
| wave 289 s, 16 s | .0697 | .0682 | .0143 (4.9×) | 19380 | .0384 | .0143 |
| wave 289 s, 64 s | .2355 | .2296 | .0843 (2.8×) | 509 | .2236 | .0843 |
| wave 5064 s, 16 s | .0220 | .0157 | .0143 (1.5×) | 18850 | .0143 | .0143 |
| wave 5064 s, 64 s | .0930 | .0843 | .0843 (1.1×) | 506 | .0843 | .0843 |
| collinear 25 Hz, 16 s | .1056 | .0860 | .0804 (1.3×) | 570 | .0804 | .0804 |
| collinear 25 Hz, 64 s | .3432 | .2939 | .2805 (1.2×) | 37.5 | .2805 | .2805 |
| sync-locked, 16 s | .1159 | .0914 | .0914 (1.3×) | 436 | .0899 | .0899 |
| sync-locked, 64 s | .3838 | .3615 | .3615 (1.1×) | 20.1 | .3615 | .3615 |

What the table shows:

- **Kill criterion.** The joint-reader `C` loses at most 4.9× on MOVING
  words, so the criterion passes; its optimum there is `k = 0`, where `C`
  does not enter.
- **Quiet words.** `kappa_W` is infinite along the kernel. The joint-reader
  loss of 43–52× at 289 s reflects the transient BA marginal, not the word;
  at 5064 s it is 5.7×.
- **Kernel.** The exact-variance kernel bound loses at most 1.9×. On steady
  quiet words the ceiling loses 1.2× with the actual tilt and 1.5× at the
  `10^-3` premise.
- **Diameter structure.** The controlling direction of `kappa_W` is tilt/BA
  (67–93% BA share), except gyro bias on the 16-s sync-locked word.
  `kappa_W` equals the slow-marginal diameter up to a fast-given-slow factor
  ≤1.078, and exceeds the persistence cap on every word.
- **Invariance.** At the `10^-3` tilt premise `kappa_nu` is 195.8–461.6
  (16 s) and 15.7–49.2 (64 s). The kernel set propagates on every word,
  checked against the next root's kernel.

These are finite replays, not enclosures.

**Reduced source-uniform obligations.** Section 6's three bounds (`J` floor,
`A - kappa J` ceiling, `I_eff` floor) reduce to:

- **(O1)** a source-uniform ceiling on the kernel-bounded diameter
  `kappa_nu`, i.e. uniform observability of every root direction except the
  one-dimensional physical tilt/BA kernel;
- **(O2)** the scalar kernel ceiling `c_nu` and its propagation
  `nu_next' P_nu nu_next <= c_next`, with the BA part already proved and a
  tilt ceiling `tau_theta` of order `10^-3 rad^2`.

No root covariance matrix, historical reader or excitation premise enters
`rho_0`. MARINE MOTION excitation only shrinks `kappa_nu`. `I_eff`/`B_*`
remain relevant to coercivity (`P <= C` in the nonlinear supplies), not to
`rho_0`.

### O1 quotient reader as an augmented minimum-action problem

Before bounding any coefficient, fix the exact inequality required by O1.
For a word (W), let the frozen design be
`y=O x_0+A_s s`, `x_N=Phi x_0+T s`, `Sigma=A_s A_s'`.
For a candidate kernel ceiling `c>0`, put `mu=1/c` and append the
fictitious scalar datum `z=sqrt(mu) nu' x_0+v`, `var(v)=1`.
Equivalently,

`O_mu=[O; sqrt(mu) nu']`, `Sigma_mu=diag(Sigma,1)`.

Its information matrix is exactly

`J_mu=J+mu nu nu'`.

Apply the joint minimum-action reader theorem to this augmented word.  The
Loewner-minimal terminal reader action is exactly

`B_mu = Pi + Phi_tilde J_mu^-1 Phi_tilde' = P_nu`.

Therefore O1 is equivalent to constructing **any** feasible quotient reader
`L_mu O_mu=Phi` whose residual action obeys

`(T-L_mu A_mu)(T-L_mu A_mu)' <= K(c) Pi`.                 (O1-R)

Indeed minimum action gives `P_nu<=B(L_mu)`, hence
`kappa_nu=lambda_max(Pi^-1 P_nu)<=K(c)`.  Conversely the minimum reader
itself attains `P_nu`, so no pivot or singular-value surrogate can improve
the mathematical target.

This identifies how the existing structure enters the reader:

- the S-chain identity first annihilates the complete `(v,p,S,a_w)` root,
  every AW process source and every AW sync in the chosen accelerometer/S
  combinations;
- applied magnetic rows read the two attitude components normal to the
  transported field direction;
- Corollary A*, G0 and Lemma T supply the remaining chronological
  attitude/gyro-bias quotient geometry;
- the appended scalar row supplies only the physical tilt/BA kernel;
- the literal injection/reset matrices stay inside `O_mu,A_mu,Phi,T`, so
  their signed temporal structure is preserved rather than norm-summed.

What remains genuinely unproved in (O1-R) is an explicit source-uniform
block reader with a finite action ratio `K(c)`: specifically the nominal
AW-window statistics and the extension of G0 to the literal injection frame.
The augmented-reader identity itself is exact and needs neither premise.

### Known-root scalar action after S-chain elimination

For O2 the exact target is not `lambda_max(Pi)` but

`d_j = nu_(j+1)' Pi_j nu_(j+1)`.

Because `Pi_j=Ric_Wj(0)`, every terminal error is a linear function of the
word sources only.  Let `q=nu_(j+1)` and write the known-root terminal
functional as

`q' x_N = q' T s`.

For any linear combination `l'y` of the same word data, the trial residual
is

`q'x_N-l'y=(q'T-l'A_s)s`

and therefore

`d_j <= (q'T-l'A_s)(q'T-l'A_s)'`.                         (D-R)

The minimum over `l` is exactly `d_j`; this is the scalar specialization
of the minimum-action reader with a known root.  Thus an explicit scalar
reader is sufficient and no full covariance upper bound is needed.

Choose the data coordinates only **after** applying the S-chain identity.
The four-S divided-difference combination cancels the complete neutral/AW
root, every AW process increment and every AW covariance sync exactly.
Consequently the residual in (D-R) contains only the attitude/BA part of the
accelerometer combinations, their measurement noises, the LIN
`(n_v,n_p,n_S)` residual already bounded in the S-chain lemma, and the
literal attitude/BA process/reset sources.  This is the correct point at
which to use the transported magnetic rows and the signed G0 geometry.
Bounding AW sources before this elimination would reintroduce DEAD_END 20.

A source-uniform number `d_bar` therefore follows once one exhibits scalar
weights `l_j` with

`sup_j |q_j'T_j-l_j'A_(s,j)|^2 <= d_bar`.                  (D-S)

The same transported quotient coordinates used in (O1-R) should be used in
(D-S); independent coefficient boxes are not admissible substitutes.  The
remaining open work is to derive those weights and bound their action from
the literal injected G0 geometry and source-uniform nominal AW-window
statistics.

### O2 reduction to a scalar known-root invariant

The separate field-axis attitude ceiling is sufficient for O2, but it is not
necessary.  The kernel-bounded diameter itself supplies a sharper route that
retains the attitude/BA correlation.

Fix a word (W_j), its root kernel `nu_j`, and a candidate scalar ceiling
`c_j>0`.  Put `mu_j=1/c_j` and
`P_nu,j = Ric_Wj(mu_j nu_j nu_j')`.  By definition,

`P_nu,j <= kappa_nu,j Pi_j`,

where `Pi_j=Ric_Wj(0)` is the known-root terminal covariance.  Corollary K
already proves that every actual root satisfying
`nu_j' P_j nu_j <= c_j` has `P_(j+1) <= P_nu,j`.  Therefore the exact
scalar propagation inequality is

`nu_(j+1)' P_(j+1) nu_(j+1)
 <= nu_(j+1)' P_nu,j nu_(j+1)
 <= kappa_nu,j d_j`,

with

`d_j := nu_(j+1)' Pi_j nu_(j+1)`.

Consequently O2 follows from the two **linked, same-word** bounds

`kappa_nu,j(1/c_j) <= K(c_j)`,  `d_j <= d_bar`,

and an invariant scalar `c_bar` satisfying

`d_bar K(c_bar) <= c_bar`.                                      (O2-FP)

This is not a full-state covariance ceiling.  `Pi_j` starts from a known
root, so `d_j` contains only the literal finite-word process/correction/reset
action.  On a bounded-duration regular A21 word its finiteness follows from
the already bounded shipping coefficients, process covariances and reset
Jacobians; what is still needed is a usable explicit source-uniform value
`d_bar`, preferably from the same block reader/S-chain construction used
for O1 rather than a norm-summed whole-word envelope.

The dependence of `K` on `c` must not be discarded.  If the physical
kernel is exactly invisible to the data, weakening the fictitious precision
`mu=1/c` lets `K(c)` grow with `c`; a bare statement
`K<infinity` does not imply (O2-FP).  The required scalar mechanism is
therefore a strict one-dimensional return inequality.  Equivalently, after
isolating the kernel coordinate in the terminal Schur complement, it is
enough to prove constants `d>=0, a>=0, b>0` such that

`c_(j+1) <= d + a c_j/(1+b c_j)`,                         (O2-R)

uniformly on recurring roots, and then exhibit `c_bar` with

`d + a c_bar/(1+b c_bar) <= c_bar`.

Here `b` is supported information in the kernel coordinate and `a` is its
literal terminal persistence.  BA decay may contribute to `a<1`; it must
not be silently replaced by an exact persistent kernel.  The proved
`P_ba,ba<=I/1600` remains available to sharpen the reader, but no independent
attitude/BA Cauchy--Schwarz split is required by (O2-FP) or (O2-R).

**Where this enters.**  (O2-FP) supplies exactly the `c` in Corollary K.
Together with O1, choose `mu=1/c_bar`,
`kappa_bar >= sup kappa_nu,j(mu)`; then the recurring linear factor is

`rho_0 <= 1 - 1/kappa_bar < 1`,  `1-rho_0 >= 1/kappa_bar`.

This margin is not yet certified because O1 and the explicit invariant
constants in (O2-FP)/(O2-R) remain open.  No nonlinear retained-radius claim
follows until those constants are substituted.

## Coupled quotient contraction and retained-radius theorem

The source-only separation between a globally uniform linear factor and a
later nonlinear radius is stronger than required.  The correct local theorem
uses one candidate radius and the same frozen actual word coefficients in
both estimates.

Let `r>0` and `c>0`.  On every regular A21 word whose every prefix satisfies
`V<=r^2` and whose root satisfies `nu'P nu<=c`, define:

- `K(c,r)`: a source-uniform upper bound on the augmented quotient diameter
  `kappa_nu`, using the radius-dependent signed AW/injection geometry;
- `D(c,r)`: a source-uniform upper bound on
  `nu_next' P_nu nu_next` (or any sharper scalar Schur return);
- `E(r)`: the complete word finite-error supply in sqrt-storage units from
  sensor/model defects, curvature, literal reset remainder, projection and
  arithmetic, with the actual-gain composition of section 3;
- `G(r)`: an every-prefix supply bound in the same units.

Corollary K gives the homogeneous word factor
`rho(c,r)<=1-1/K(c,r)`.  Put `q(c,r)=sqrt(rho(c,r))`.
The exact finite-error composition then gives

`sqrt(V_next) <= q(c,r) sqrt(V_root) + E(r)`.              (CR-1)

The kernel set propagates simultaneously when

`D(c,r) <= c`.                                             (CR-2)

Thus the rectangle
`R(r,c)={V<=r^2, nu'P nu<=c}` is invariant at recurring roots provided

`D(c,r)<=c,`
`q(c,r) r + E(r) <= r`.                                   (CR-3)

Every-prefix retention is obtained from the already exact prefix composition
when

`r + G(r) <= r_prefix`

for a declared prefix radius `r_prefix` on which all coefficient, projection
and reset bounds used in `K,D,E` remain valid, and the next complete word
returns to `R(r,c)`.  In the convenient no-enlargement form this is simply
a direct prefix self-map inequality for each prefix.

Equivalently the strict complete-word margin is

`M(r,c) := (1-q(c,r))r-E(r) > 0`,                          (CR-M)

together with the scalar-kernel margin `c-D(c,r)>=0`.
This is the quantity that must be positive before any retained radius is
claimed.

There is no circularity: choose candidate `(r,c)`; derive `K,D,E,G`
uniformly over that closed candidate set from the literal shipping
operations; then verify (CR-3) and the prefix inequalities.  The resulting
invariance proves that the assumptions used to derive the bounds remain true
on subsequent words.

The radius-dependent geometry has the schematic form

`m_perp(r)<=m0+m1 r+m2 r^2`,
`delta_Q(r)<=q0+q1 r+q2 r^2+q3 r^3`.

These enter the local-tube quotient reader directly, producing `K(c,r)` and
`D(c,r)`; they are not promoted to source-only constants.  The accelerometer
curvature `|f_hat||theta|^2/2+|theta||e_aw|`, magnetic curvature, finite
reset remainder and projection/arithmetic terms enter `E(r)` exactly where
section 3 places their operation defects.

If (CR-3) closes, the certified local linear margin on the invariant set is

`1-rho_0(r,c) >= 1/K(c,r) > 0`,

but regional practical stability still additionally requires startup/H18
entry, regime transitions, recurring magnetic-service qualification and
implementation/float32 totality.

## Canonical analytic local-tube readers and actions

Let a local-tube/S-chain reduced record be
`z=H x + A_z s`, with root quotient coordinate `x=(theta,beta)`, source
factor `s`, and positive residual covariance `Sigma_z=A_z A_z'`.
Append the fictitious kernel row with precision `mu=1/c`.  Put

`H_mu=[H; sqrt(mu) nu']`, `Sigma_mu=diag(Sigma_z,1)`,
`G_mu=H_mu' Sigma_mu^-1 H_mu`.

For any terminal quotient map `F_q`, the canonical weighted reader is

`L_q = F_q G_mu^-1 H_mu' Sigma_mu^-1`.                    (LR-1)

Its coefficients are therefore explicit functions of the literal transported
tube rows.  It satisfies `L_q H_mu=F_q`.  Its residual source coefficient
and quadratic action are

`Z_q=T_q-L_q A_mu`,
`R_q=Z_q Z_q'`.                                            (LR-2)

The minimum-action reader can only improve this action, so (LR-2) is a valid
constructive O1 reader whenever `G_mu>0`.

For the scalar next-kernel functional `q'=nu_next'`, use

`l_d = q' T_x G_mu^-1 H_mu' Sigma_mu^-1`,                  (LR-3)

where `T_x` is the root-to-terminal quotient map.  Its source residual and
action are

`z_d=q'T_s-l_d A_mu`,
`R_d=z_d z_d'`.                                             (LR-4)

Thus `D(c,r)<=R_d` directly.  Equations (LR-1)--(LR-4) are the requested
analytic reader coefficients and quadratic actions; they retain all source
cross terms through the common factor `A_mu`.

The local G0 inequalities give an explicit lower matrix for `G_mu`.
With `a=K_II gamma`, `b=K_II sqrt(gamma) C_II`, and the independent gyro
supply `q_I`, the radial two-block comparison is

`G_mu >= [[a I,-b I],[-b I,(q_I+K_II C_II^2)I]] + mu nu nu'`.   (LR-5)

Equivalently, before the kernel row, the Schur complement in the gyro block is
at least `q_I`; the combined smallest quotient eigenvalue is the G0
`s(c,r)^2`.  Hence

`||G_mu^-1|| <= 1/s(c,r)^2`                                (LR-6)

on the quotient.  This yields the coarse explicit action estimate

`R_q <= (||T_q|| + ||F_q|| ||A_mu||/(s(c,r) sqrt(lambda_min(Sigma_mu))))^2 I`, (LR-7)

and the analogous scalar estimate

`R_d <= (||q'T_s|| + ||q'T_x|| ||A_mu||/(s(c,r) sqrt(lambda_min(Sigma_mu))))^2`. (LR-8)

These norm ceilings are existence bounds only and may reproduce DEAD_END 22
if used globally.  The useful action is the exact signed quadratic form
(LR-2)/(LR-4), evaluated symbolically after the S-chain cancellations and
tube-to-tube coefficient transport, not (LR-7)/(LR-8).

The missing transport constants now have exact definitions rather than free
parameters.  For the normalized force weights `alpha_k` on an anchored tube,

`A0=|sum alpha_k u_phys,k x b|`,
`A1=sup_{V<=r^2}|sum alpha_k (u_hat,k-u_phys,k) x b|/r`,

with the error term expanded by the literal AW-loop identity before taking
the supremum.  For relative reset transport `Q_j Q_k'`, write each literal
factor as `I-X_l/2+R_l`, `||R_l||<=|x_l|^3/6`.  Then `Q0` is the
integrated physical gyro-residual/implementation term, `Q1` the coefficient
of the signed partial sum `S_n=sum_l x_l`, `Q2` the exact pair-product
coefficient `sum_l |x_l||S_(l-1)|/4`, and `Q3` the accumulated
`sum_l |x_l|^3/6` remainder after the same signed transport.  No independent
per-reset box is introduced.

What remains unproved is not the reader formula: it is a finite
source/radius-uniform evaluation of these signed suprema and of the exact
actions (LR-2),(LR-4).  Until those are bounded analytically, replacing them
by carried-word values would be fitted rather than rigorous.

## Signed factor-space bounds for R_q and R_d

The actions (LR-2),(LR-4) admit a nonmultiplicative analytic bound after the
S-chain reduction.  Index the independent fresh source blocks by `a` and
write the reduced record and terminal maps as

`A_mu=[A_a]_a`, `T_q=[T_qa]_a`, `T_s=[T_sa]_a`.

The same source block may enter many rows; it is represented once.  For the
canonical reader define the signed factor coefficients

`Z_qa=T_qa-L_q A_a`,
`z_da=q'T_sa-l_d A_a`.                                      (FA-1)

Because the whitened fresh factors are mutually independent by construction,

`R_q=sum_a Z_qa Z_qa'`,
`R_d=sum_a z_da z_da'`.                                     (FA-2)

No cross-factor absolute-value bound is needed.  Correlations within one
prediction factor, including attitude/gyro-bias and LIN/AW cross terms, stay
inside its matrix block `A_a,T_a` and are squared only after the signed
cancellation in (FA-1).

For an S-chain spanning steps k, the AW/root/sync columns are exactly zero.
Its remaining LIN factor coefficient is
`-g_k'(n_v,n_p,n_S)_k`; hence its contribution to the reduced-record
covariance is bounded by

`V_S <= (1+eps) sum_k q_k |g_k|^2
               (h_k^3/3+h_k^5/20+h_k^7/252)`,               (FA-3)

with `q_k<=2 sigma_k^2/tau_k`.  Measurement factors contribute their actual
bounded `R_acc,R_mag,R_S` after the same signed row weights.  BA and AG
prediction factors remain chronological blocks and are not split into
independent scalar boxes.

Now partition the word into anchored tubes r.  Let `C_r` be the exact
invertible quotient coordinate map from tube r to a common terminal anchor.
Transport the local reader coefficient, not its norm:

`L_(q,r)^term=C_r L_(q,r)`,
`l_(d,r)^term=l_(d,r) C_r^-1`.                               (FA-4)

For each physical source block a, sum all of its transported occurrences
before squaring:

`Xi_qa=T_qa-sum_r L_(q,r)^term A_(r,a)`,
`xi_da=q'T_sa-sum_r l_(d,r)^term A_(r,a)`.                   (FA-5)

Then the complete-word actions are exactly

`R_q=sum_a Xi_qa Xi_qa'`,
`R_d=sum_a xi_da xi_da'`.                                   (FA-6)

This is the desired signed tube-to-tube action formula.  It cannot exhibit
the exponential operation-count growth of the retired pivot bound: the
number of terms is the number of fresh source blocks, and cancellation across
all rows/tubes sharing a block occurs in `Xi` before its norm is taken.

A rigorous scalar ceiling follows from per-family energy bounds without
destroying temporal signs.  If source blocks are grouped into disjoint
families `F` (AG process, BA process, S-chain LIN, acc noise, magnetic
noise, S noise, implementation defect), define

`B_(q,F)=sum_(a in F) Xi_qa Xi_qa'`,
`b_(d,F)=sum_(a in F) xi_da xi_da'`.                         (FA-7)

Then

`R_q=sum_F B_(q,F)`, `R_d=sum_F b_(d,F)`,                 (FA-8)

and it is sufficient to prove matrix/scalar ceilings
`B_(q,F)<=bar B_(q,F)`, `b_(d,F)<=bar b_(d,F)`.
The already proved S-chain bound (FA-3) closes the LIN family.  The BA family
has total fresh variance bounded by its literal OU/RW recursion and the
proved marginal `P_ba<=I/1600`.  Measurement families have fixed shipping
noise ceilings.  The remaining unclosed family is AG process plus its
signed injected transport: its `Xi` contains the same local-tube
`Q_k,D_k` coefficients that determine G0.  Thus the action problem has been
reduced to one matrix-energy bound on that family; AW/root/sync and nuisance
families no longer obstruct it.

Consequently define

`R_q^bar(c,r)=sum_F bar B_(q,F)(c,r)`,
`R_d^bar(c,r)=sum_F bar b_(d,F)(c,r)`.                       (FA-9)

Then the rigorous coupled coefficients are

`K(c,r)=lambda_max(Pi_lower^-1 R_q^bar(c,r))`,
`D(c,r)=R_d^bar(c,r)`.                                      (FA-10)

Here `Pi_lower` is the already proved regular-root lower covariance, with
the same coordinates.  The only new bound still required for numerical
closure is the AG-family signed energy in (FA-7); replacing it by
`||C_r||||L_r||||A_a||` separately would destroy the cancellation and
recreate DEAD_END 22.

## AG-process signed energy: inverse-frame inequality

The remaining family in (FA-7) can be bounded without forward products of
reset norms.  Work in the local anchor inverse frame.  Let the chronological
AG map be `T_k=[[A_k,B_k],[0,I]]` and put `C_k=A_k^-1 B_k`.  A literal
reset `G=I+[d]/2` sends `A,B` to `GA,GB`, hence leaves `C` unchanged
exactly.  Moreover `sigma_min(G)>=1`, so `||A_k^-1||<=1` for a word
started at an anchor with `A=I`; Rodrigues predictions are orthogonal and
the qualified small-rate branch has the same inverse nonexpansion already
proved in the chronological-transport lemma.

For prediction k write the full correlated AG process factor as
`U_k=[U_theta,k;U_b,k]`, `Q_k=U_k U_k'`.  In inverse-frame quotient
coordinates its fresh contribution is

`W_k=[[A_(k+1)^-1, -C_(k+1)],[0,I]] U_k`.                  (AG-1)

The minus sign is essential: attitude and gyro-bias process columns are not
split.  The complete signed reader coefficient for this fresh block is

`Xi_(q,k)=H_(q,k) W_k`, `xi_(d,k)=h_(d,k) W_k`,           (AG-2)

where `H_(q,k)`, `h_(d,k)` are the already transported residual-reader
maps from the injection time to the terminal anchor after all row
cancellations.  Therefore

`B_(q,AG)=sum_k H_(q,k) W_k W_k' H_(q,k)'`,
`b_(d,AG)=sum_k h_(d,k) W_k W_k' h_(d,k)'`.                (AG-3)

This is exact factor-space energy.

A source-uniform Loewner bound follows directly from the literal Q blocks.
Let
`Q_tt,k<=q_theta,k I`, `Q_bb,k=q_b h_k I`, and
`||Q_tb,k||<=q_x,k`.  Since `||A^-1||<=1`,

`W_k W_k' <= E_k(C_(k+1))`,                                (AG-4)

with the explicit 6x6 block majorant

`E_k(C) =
 [[q_theta,k I + q_b h_k C C' + q_x,k( C+C')_abs,  *],
  [*, q_b h_k I]]`.

For a scalar safe form, completing the square with any `eta_k>0` gives

`W_k W_k' <= diag(
 (1+eta_k) q_theta,k I
 +(1+1/eta_k) q_b h_k C C'
 +(1+eta_k) 2 q_x,k I,
 (q_b h_k+2 q_x,k/eta_k) I )`.                              (AG-5)

The useful bound keeps the matrix `C C'` and optimizes `eta_k` only after
the signed reader coefficient `H_(q,k)` is applied.

The inverse-frame recurrence supplies

`C_(k+1)=C_k+A_k^-1 R_k^-1 D_k` at predictions and no change at resets.
Hence for any anchor interval

`C_k=C_0+sum_(j<k) A_j^-1 R_j^-1 D_j`.                     (AG-6)

Every summand has norm at most `h_j(1+e_D)` on the qualified prediction
domain.  More importantly, (AG-6) is a signed chronological sum; it is the
same `D_k` history used by the local-tube G0 reader.  Substitute (AG-6)
into (AG-3) before taking norms.  Define the tube energy matrices

`M_(q,theta)=sum_k H_(q,k)H_(q,k)'`,
`M_(q,C)=sum_k H_(q,k) C_(k+1)C_(k+1)' H_(q,k)'`,          (AG-7)

and analogously `m_(d,theta),m_(d,C)`.  Then (AG-5) gives

`B_(q,AG) <= a_theta M_(q,theta)+a_C M_(q,C)+a_b M_(q,b)`, (AG-8)
`b_(d,AG) <= a_theta m_(d,theta)+a_C m_(d,C)+a_b m_(d,b)`, (AG-9)

where the coefficients are explicit sums/maxima of the shipping
`q_theta,k,q_b h_k,q_x,k` and the chosen Young parameters.  No reset-product
factor appears.

Equations (AG-6)--(AG-9) are the controlling matrix-energy inequality.  The
only remaining quantity is the signed reader-weighted chronological energy
`M_(q,C)` (and its scalar d analogue).  It is not an independent new
premise: `C_k` is exactly the gyro-history coordinate already constrained
by the local-tube G0/Lemma-T construction.  A tube bound
`sum_k H_k C_k C_k' H_k'<=C_E(c,r)` therefore closes the last source family
and yields finite `R_q^bar,R_d^bar`.

This result removes literal reset amplification and preserves the correlated
AG process factor.  It does not yet assign a numerical `C_E`; replacing
the signed energy by `sum ||H_k||^2 ||C_k||^2` over a long word is the
coarse norm route and is not promoted.

## Reader-weighted C-energy: leverage inequality

Lemma T supplies lower information, not an upper bound on `|C_k|`; therefore
it must enter through the inverse normal matrix of the same tube reader.
Let the whitened reduced tube rows be `X_k` and stack `X=[X_k]_k`.
Include the fictitious kernel row and put

`G_mu=X'X+mu nu nu'`.

For terminal quotient map `F`, the canonical reader block on row k is

`L_k=F G_mu^-1 X_k'`.                                      (LE-1)

Let `C_k` denote the inverse-frame gyro-history coordinate at that row and
define the block-diagonal multiplication operator
`mathcal C=diag(C_k)`.  The signed reader-weighted history energy is exactly

`M_C=F G_mu^-1 X' mathcal C mathcal C' X G_mu^-1 F'`.      (LE-2)

Do not bound `mathcal C` separately.  Decompose the tube rows into attitude
and gyro columns, `X_k=[U_k,U_k C_k+V_k]`, where `V_k` contains the
literal within-step gyro injection and the transported magnetic contribution.
Then

`X [0;I] = mathcal U C + mathcal V`

in stacked notation.  The G0/Lemma-T proof gives a positive lower quadratic
form on this same gyro column after eliminating attitude.  Denote its Schur
complement by

`S_g = X_g'(I-P_U)X_g + mu S_nu >= q_T I`,                 (LE-3)

where `P_U` is the weighted attitude-column projector and `q_T>0` is the
literal local-tube Lemma-T/G0 gyro floor (including the kernel row when
needed).

The weighted least-squares leverage inequality applies only to the residualized gyro columns

`Xg_perp=(I-P_U)X_g`,
`Xg_perp S_g^-1 Xg_perp' <= I`                              (LE-4 corrected)

(up to the kernel-row augmentation in the same residualized space). The previously stated unprojected inequality `X_g S_g^-1 X_g'<=I` is false. Consequently the following AG-process action must be rederived with `Xg_perp` carried consistently. Writing the AG-process coefficient as observed gyro-column part plus the within-step remainder,

`mathcal C = mathcal C_obs + mathcal R_D`,

the canonical reader action obeys

`M_C <= ...`                                                (LE-5 RETRACTED)

The former displayed bound is not established by LE-4 corrected because its decomposition used the unprojected gyro column.

The first term is controlled directly by the quotient information:
`F G_mu^-1 F'<=||F||^2/s(c,r)^2 I`.  The second term contains only the
within-step defect `R_D`, not the accumulated chronological `C_k`.
The implemented gyro invariant already gives
`||R_k^-1D_k-h_k I||<=h_k(theta_k/2+theta_k^2/3)`; hence with
`theta_k<=theta_max<.007`,

`||R_D,k||<=h_k e_D`,
`e_D=theta_max/2+theta_max^2/3`.                            (LE-6)

Therefore LE-7 and LE-8 are RETRACTED. A valid AG-process energy bound requires re-expressing the source coefficient in the residualized attitude-orthogonal gyro coordinates and separately charging the `P_U X_g` component. No later theorem may cite LE-5--LE-8 as established.

Thus the accumulated signed `C_k` energy is bounded without
`sum ||H_k||^2||C_k||^2` and without a reset-product norm.  Lemma T enters
only through the positive Schur floor `q_T`; large chronological gyro
history increases both the raw coefficient and the information that
normalizes its reader leverage.

Scope: (LE-5)--(LE-8) require the algebraic decomposition of the literal
S-chain/local-tube rows into the same `X_g` used by the G0 Schur complement.
The injection/reset frame must therefore use one common whitening and anchor
convention.  This identification is exact for the ideal local-tube array and
remains to be checked for the literal finite-series/implementation defect.
Until that bookkeeping check is closed, (LE-7) is a conditional analytic
bound, not a numerical shipping certificate.

## Literal whitening and anchor identification for the leverage bound

The bookkeeping identification required by (LE-3)--(LE-8) is exact in real
arithmetic provided whitening is performed **after** the S-chain/local-tube
row combinations.

Shipping conventions are:

1. the mean quaternion is WORLD-to-BODY and prediction uses the normalized
   `quat_from_delta_theta(-w h)`;
2. the AG covariance prediction is
   `[[R_s,B_s],[0,I]]` from the same `rot_and_B_from_wt(w,h)`;
3. an applied correction changes the mean by the literal quaternion injection
   and applies covariance reset `G=I+[d]x/2`;
4. the world proof coordinates are obtained by left multiplying body rows by
   the orthogonal `R_k'`, giving
   `[f_k]x[Atilde_k,Btilde_k]` exactly as in Lemma W.

Thus the local anchor convention used in the leverage proof is

`Q_k=Atilde_k Atilde_j^-1`,
`D_k=Atilde_j^-1(Btilde_k-Btilde_j)`,

with the same chronological reset factors as shipping.  No transpose or
sign change is missing: `rot_and_B_from_wt` is the WORLD-to-BODY covariance
transition, while the proof's `R_k'` is precisely the body-to-world
orthogonal change of row coordinates.

Finite-series prediction terms do not alter this identification.  They are
already inside the literal `R_s,B_s`.  Relative to the ideal integral they
appear only through

`R_s^-1 B_s = h I + Delta_D`,
`||Delta_D||<=h(theta/2+theta^2/3)`

on the qualified small-angle domain.  This is exactly the `R_D` defect in
(LE-5)--(LE-7), not a change of anchor or whitening.

The whitening point is essential.  Let `C_S` denote the deterministic
matrix that forms all selected S-chain and local-tube combinations from the
raw observation vector.  If the raw auxiliary record is

`y=O x+A s`, `Sigma=A A'`,

then the reduced record is

`z=C_S y=H x+A_z s`,
`H=C_S O`, `A_z=C_S A`,
`Sigma_z=C_S Sigma C_S'=A_z A_z'`.                         (WH-1)

The S-chain uses common prediction factors in several combined rows, so
`Sigma_z` is generally **not block diagonal** even though the fresh factors
are independent.  Therefore the literal whitened array for Lemma T/G0 is

`X=Sigma_z^-1/2 H`,                                        (WH-2)

with any common square root/inverse factor of the full reduced covariance.
Whitening individual accelerometer, magnetic or S rows before applying
`C_S` is not equivalent and is not allowed.

Under (WH-2),

`X'X=H' Sigma_z^-1 H`

is exactly the reduced-record information used by the canonical reader.
Partitioning its columns into attitude and gyro coordinates therefore gives
the same `X_g` and the same weighted attitude projector `P_U` that appear
in the Schur complement (LE-3).  The fictitious kernel row is appended after
this reduction with independent variance one, so it adds exactly
`mu nu nu'`.

Consequently the **residualized** leverage identity

`Xg_perp S_g^-1 Xg_perp'<=I`,  `Xg_perp=(I-P_U)X_g`,

uses the same literal whitening and anchor convention as the reader action.
The unprojected statement `X_g S_g^-1 X_g'<=I` is false, as already noted
above, and is not restored by this bookkeeping identification.  The literal
reset transport is retained exactly; finite-series prediction error is
isolated in `R_D`; S-chain shared-source correlations are retained in
`Sigma_z`.

**Result.**  The whitening/anchor identification passes only for the
residualized operator.  It does **not** establish the former signed
`C_k`-energy bound.  The missing `P_U X_g` component must be charged
separately, so LE-5--LE-8 and the claimed AG-process family bound remain
RETRACTED.

## Direct literal reduced-information formulation

The AW-mean and reset-product scalar pre-bounds can be removed entirely.
Start from the raw frozen auxiliary record

`y=O x+A s`, `Sigma=A A'`.

Choose only deterministic S-chain/local-tube combinations `C_S` whose
coefficients depend on timestamps and literal transition coefficients, not on
unknown source values.  Form

`H=C_S O`, `A_z=C_S A`, `Sigma_z=A_z A_z'`.             (DI-1)

All AW root columns, AW process columns and AW covariance-sync columns selected
by the S-chain vanish **algebraically** in `H,A_z`; no nominal AW mean bound
is used.  Accelerometer attitude columns remain with their actual nominal
`a_hat_k` values as linked coefficients of the frozen word.

Retain every literal reset factor in the historical AG columns.  In world
coordinates these are the exact factors

`N_l=R_(l+) ' G_l R_(l-)`,
`N_l'N_l=I+(|x_l|^2 I-x_l x_l')/4 >= I`.                  (DI-2)

Do not replace their ordered product by a rotation plus Q1/Q2/Q3 remainder.
The inequality `sigma_min(N_l)>=1` is useful for invertibility and inverse
nonexpansion, but by itself does not preserve force/field kernel angles;
therefore no scalar Corollary-A transport through the reset product is
claimed.

Whiten only after reduction:

`X=Sigma_z^-1/2 H`.                                        (DI-3)

Partition the slow columns as physical kernel coordinate `nu` and a chosen
quotient basis `E_q`.  Append the fictitious kernel row and define the exact
reduced normal matrix

`G_red(c)=
 [E_q,nu]' H' Sigma_z^-1 H [E_q,nu]
 +diag(0,1/c)`.                                             (DI-4)

Equivalently avoid a basis by using the projector onto `nu^perp`.  The
literal quotient and gyro Schur floors are

`s_lit(c)^2=lambda_min(G_red(c))`,
`q_T,lit(c)=lambda_min(Schur_gyro(G_red(c)))`.              (DI-5)

Every `N_l`, actual nominal AW coefficient, finite-series `R_s,B_s`,
measurement covariance and shared process correlation is inside (DI-4).
Thus A0/A1 and Q0..Q3 disappear from the theorem.

For a radius-local proof, coefficients generated by estimator errors range
over the candidate set `R(r,c)`.  The actual obligation is the matrix
enclosure

`G_red(c;history) >= G_*(c,r)>0
 for every history in R(r,c)`.                              (DI-6)

This is a structured same-history matrix inequality, not independent boxes.
MARINE MOTION and MAGNETIC SERVICE enter only to prove (DI-6), for example by
pairing actual reduced accelerometer rows with service magnetic rows and the
inverse-frame chronological gyro coordinate.  No separate nominal-AW mean or
reset-product estimate is required.

Once (DI-6) is proved, set

`s(c,r)^2=lambda_min(G_*(c,r))`,
`q_T(c,r)=lambda_min(Schur_gyro(G_*(c,r)))`,

and all leverage/factor-space action bounds above apply unchanged.

**Consequence.**  The previous A_i/Q_i route is retired as unnecessary.
The controlling unresolved statement is now one direct lower Loewner bound
(DI-6) for the literal correlated reduced information matrix.  It is not
proved merely by `sigma_min(N_l)>=1`: noncontractive resets can still change
kernel alignment.  A proof of (DI-6) must exploit the complete linked row
matrix and physical/service chronology.  This is narrower and avoids the
circular cumulative correction-action premise.

## Variational reduced-information floor

The direct lower bound should use the **optimal nuisance-annihilating
compression**, not an arbitrary S-chain row selection.

Partition the raw frozen design into slow/root columns `O_s` and canceled
fast/nuisance columns `O_f`:

`y=O_s x_s+O_f x_f+A s`, `Sigma=A A'>0`.

For fixed `x_s`, minimize the whitened data energy over the unrestricted
nuisance root `x_f`:

`Q_red(x_s)=min_(x_f) ||Sigma^-1/2(O_s x_s+O_f x_f)||^2`.   (VI-1)

The minimizer is the weighted projection onto the complement of
`range(Sigma^-1/2 O_f)`, hence

`G_red =
 O_s' Sigma^-1 O_s
 -O_s' Sigma^-1 O_f
  (O_f' Sigma^-1 O_f)^dagger
  O_f' Sigma^-1 O_s`.                                      (VI-2)

Equivalently, if `P_f` is the orthogonal projector onto
`range(Sigma^-1/2 O_f)`,

`G_red=O_s' Sigma^-1/2 (I-P_f) Sigma^-1/2 O_s`.            (VI-3)

This is the maximal information obtainable by any linear compression that
annihilates `O_f`.  Every explicit S-chain reader is a feasible compression
and therefore lies below (VI-2); it need not preserve a source-uniform floor.
The S-chain identities remain useful to prove that the canceled
`(v,p,S,a_w)` root/AW sources belong to the nuisance span, but the theorem
should use the Schur complement (VI-2) itself.

Append the physical-kernel precision after nuisance elimination:

`G_red,mu=G_red+mu nu nu'`, `mu=1/c`.                    (VI-4)

The controlling obligation is therefore the variational inequality

`Q_red(x_s)+mu(nu'x_s)^2 >= g_*(c,r)|x_s|^2`              (VI-5)

for every slow vector and every admissible same-history word in the candidate
region.  Then `G_*(c,r)=g_* I` is a valid DI-6 floor; block versions may be
used to retain a sharper gyro Schur constant.

The physical chronology enters (VI-5) without independent coefficient boxes:

- applied magnetic rows penalize the two attitude components transverse to
  the actual transported field;
- accelerometer rows and their exact nuisance projection penalize the
  complementary attitude/BA combinations;
- MAGNETIC SERVICE bounds gaps between actual magnetic rows;
- chronological gyro transport/Lemma T prevents a nonzero gyro-bias quotient
  from remaining in the moving field-axis null line across separated service
  tubes;
- MARINE MOTION attitude-span supplies the physical change needed on every
  complete moving excitation window; STILL is handled by its separate
  compatible-class argument.

A contradiction proof of (VI-5) is now natural.  If no positive uniform floor
exists on a compact candidate history class, take a sequence of unit slow
vectors with `Q_red+mu(nu'x)^2->0`.  The residuals imply, successively:
magnetic rows force attitude into the transported field-axis line; nuisance-
projected accelerometer rows force the associated BA/tilt compatibility;
Lemma T plus recurring service forces gyro bias into the compatible axial
component; the fictitious kernel row then removes the sole remaining physical
tilt/BA line.  The limit slow vector must be zero, contradicting unit norm.

What remains for a **quantitative** floor is to turn those four implications
into explicit inequalities with constants on the radius-local compact class.
The exact reset factors and literal nominal AW values remain inside
`O_s,O_f,Sigma`; no A_i or Q_i bounds are required.  Compactness alone would
prove existence of `g_*>0`, but a usable `K,D,E` requires explicit
moduli for the magnetic, nuisance-projected accelerometer and Lemma-T steps.

## Closed-loop null invariance and radius-local reachability target

Actual optimal gains cannot repair an exact raw Jacobian null direction.  At
any correction,

`e_plus=(I-KH)e_minus`.

Hence `H v=0` implies `(I-KH)v=v` for every realized gain `K`; Joseph
conditioning changes covariance but not this deterministic homogeneous
direction.  Interleaved gains can help only after prediction/reset transport
moves the direction out of the next row nullspace.  Therefore the missing
physical-to-nominal bridge cannot be obtained from gain magnitude or NIS
alone.

The remaining viable mechanism is radius-local reachability.  Let the true
world inertial acceleration be `a(t)`, true residual accelerometer bias
`b_a(t)`, nominal states `a_hat(t),b_hat_a(t)`, and attitude error
`R_hat=R_err R_true`.  On `V<=r^2`, covariance coercivity gives

`|delta theta|<=c_theta(r) r`,
`|a_hat-a|<=c_aw(r) r`,
`|b_hat_a-b_a|<=c_ba(r) r`,                                (RE-1)

**only if** the corresponding covariance upper bounds are available on the
candidate set.  The AW marginal ceiling gives the one-sided error estimate
`|e_aw|<=4r`; the BA marginal gives `|e_ba|<=r/40`.  Thus the nominal
world force differs from the true world force by at most `4r` in the AW
coordinate, without invoking pointwise AW tracking as a physical premise:
it is a consequence of being inside the candidate storage ball.

For any unit field direction `b` and any normalized convex window weights,

`| (sum alpha a_hat) x b |
 <= |(sum alpha a) x b| + 4r`.                              (RE-2)

More importantly for a lower separation, if the physical window supplies

`| (sum alpha a) x b | >= m_phys`,

then

`| (sum alpha a_hat) x b | >= m_phys-4r`.                  (RE-3)

This is a radius-local reachability exclusion of nominal collinearity whenever
`r<m_phys/4`.  It uses the proved AW covariance marginal and the definition
of storage, not an independent source assumption.

However MARINE MOTION's current attitude-span premise does **not** imply a
positive `m_phys` for the signed acceleration mean.  A vessel may change
attitude while its translational inertial acceleration has zero or
field-parallel weighted mean.  Therefore (RE-3) alone does not close the
moving information floor.

The correct reachability target must use gravity direction in the body-frame
measurement rather than inertial-acceleration mean.  For two times whose true
attitudes differ by at least `Delta_R`, the body gravity directions differ
by at least `2 sin(Delta_R/2)`.  Candidate attitude error perturbs each by at
most `2 sin(c_theta(r)r/2)`.  Thus the nominal body gravity directions retain
separation

`Delta_g,nom(r) >=
 2 sin(Delta_R/2)-4 sin(c_theta(r)r/2)`.                    (RE-4)

If positive, two literal accelerometer attitude Jacobians cannot share the
same gravity-axis kernel after nuisance projection **unless** AW/BA nuisance
columns mimic their difference.  Bounding that mimic is now a finite
two-epoch Schur problem using the proved AW/BA covariance/action bounds, not a
pointwise AW tracking statement.

So the new quantitative target is the two-epoch nuisance Schur inequality

`Q_acc,red(theta,ba)
 >= a_2(r) dist((theta,ba),K_phys)^2`, `a_2(r)>0`,         (RE-5)

with `a_2(r)` derived from (RE-4), Racc, the AW/BA nuisance action, and the
literal two-epoch transition.  Magnetic service then intersects `K_phys`
with the field-axis line, Lemma T handles gyro bias, and the fictitious kernel
row removes the remaining one-dimensional compatibility class.

This route uses the existing MARINE MOTION attitude-span assumption directly
and avoids requiring a physical acceleration mean.

## Literal two-epoch accelerometer/AW/BA Schur modulus

At two applied accelerometer epochs, rotate residuals into common world
coordinates and write the slow attitude/BA contribution as
\`d(z)=D z\`, \`z=(theta,b_a)\`, with

\`D=[[C_0,B_0],[C_1 T_theta,B_1 phi_b]]\`.

Here \`C_i=-[f_i]x\` are the literal nominal specific-force attitude rows,
\`B_i\` are the transported BA rows, and \`T_theta\` is the literal
inter-epoch attitude transport.  The AW nuisance obeys
\`a_1=Phi_a a_0+w_a\`.  If \`U_a\` bounds the inherited AW-root action and
\`Q_a\` the inter-epoch AW process action, define

\`N=[[I,0],[Phi_a,I]]\`,
\`S_a=N diag(U_a,Q_a) N'+diag(R_0,R_1)\`.

Eliminating the AW root and process factor gives exactly

\`Q_acc,red(z)=d(z)' S_a^-1 d(z)\`,
\`G_2=D' S_a^-1 D\`.

Therefore, with \`K_phys=ker D\`,

\`Q_acc,red(z)>=a_2 dist(z,K_phys)^2\`,
\`a_2=lambda_min^+(D' S_a^-1 D)\`.

A source-valid explicit lower bound is

\`a_2 >= sigma_min^+(D)^2/lambda_max(S_a)\`,

and

\`lambda_max(S_a)
 <= R_acc,max+(1+|Phi_a|)^2 U_a+Q_a,max\`.

Thus

\`a_2(r) >= d_2(r)^2/
 [R_acc,max+(1+|Phi_a|)^2 U_a+Q_a,max]\`

whenever a source-uniform geometric bound
\`sigma_min^+(D)>=d_2(r)>0\` is proved.

This closes the AW/BA Schur algebra: unrestricted AW nuisance does not erase
information for a fixed nondegenerate D; it only enlarges the finite
denominator through its root/process action.

The remaining issue is the numerator.  MARINE MOTION attitude span does not
by itself imply \`d_2(r)>0\`: the literal \`C_i\` contain nominal specific
force, not gravity alone, and admissible translational acceleration can
compensate a change in gravity direction.  The radius estimate
\`|a_hat-a|<=4r\` transfers a physical specific-force separation if one is
available, but the current attitude-span assumption does not supply that
separation.  Hence no positive numerical \`a_2(r)\` is claimed from attitude
span alone.

## Two-epoch geometric numerator: exact elimination and obstruction

For
\`D=[[C_0,B_0],[C_1 T_theta,B_1 phi_b]]\`
with orthogonal transported BA rows \`B_i\`, eliminate \`b_a\` from the
noise-free two-epoch equations.  The first row gives
\`b_a=-B_0' C_0 theta\`.  Substitution into the second gives the exact relative
operator

\`L_2=C_1 T_theta-phi_b B_1 B_0' C_0\`.                      (GE-1)

Hence
\`ker D\` is isomorphic to \`ker L_2\`, and the nonzero singular geometry of D
is controlled by L_2.  A quantitative comparison follows from the elimination
map:
for any z=(theta,b_a), let \`e_0=C_0 theta+B_0 b_a\`.  Then

\`C_1T_theta theta+B_1 phi_b b_a
 =L_2 theta+B_1 phi_b B_0' e_0\`.

Thus
\`|Dz|^2=|e_0|^2+|L_2 theta+B_1 phi_b B_0'e_0|^2\`.
Completing the square yields, for any eta>0,

\`|Dz|^2 >= [eta/(1+eta)] |L_2 theta|^2
            +[1-eta phi_b^2] |e_0|^2\`,                     (GE-2)

whenever \`eta<1/phi_b^2\`.  Therefore an explicit d2 follows from a positive
singular floor of L2 together with the bounded reconstruction
\`b_a=B_0'(e_0-C_0 theta)\`.

But current MARINE MOTION attitude span does not imply
\`sigma_min^+(L_2)>0\`.  The matrices \`C_i=-[f_i]x\` depend on nominal
specific force.  It is admissible for translational acceleration to make the
two transported nominal force directions satisfy

\`C_1 T_theta=phi_b B_1B_0' C_0\`

on the relevant two epochs, in which case \`L_2=0\` and D has only the
single-epoch rank.  Physical attitude can change between the epochs while the
specific-force cross-product maps remain aligned because translational
acceleration compensates gravity.

The existing same-cell/collinear constructions in the world-frame proof show
this mechanism explicitly: bounded displacement, velocity, acceleration,
jerk and nonzero attitude span can coexist with force parallel to the magnetic
direction at selected correction times.  MAGNETIC SERVICE excludes one
particular integer-second cadence but does not convert attitude span into a
uniform two-epoch specific-force separation.

Consequently no source-uniform \`d_2(r)>0\` can presently be derived from the
existing MARINE MOTION attitude-span assumption plus reachability radius
alone.  The two-epoch route is therefore CLOSED as a controlling path unless
another already-proved shipping-history condition supplies relative
specific-force excitation.  Do not continue trying to obtain the A21 floor
from attitude span alone.

This does not refute the full recurring A21 theorem: magnetic-service rows,
more than two accelerometer epochs, S pseudo-observations and closed-loop
chronology may jointly remove the degeneracy.  It refutes only the proposed
two-epoch accelerometer bridge as a source-uniform consequence of the current
physical assumptions.

## Complete multi-epoch corrected-word operator

The failed two-epoch reduction should not be iterated.  Use the complete
regular A21 word as one conditional Gaussian operator.

Freeze the literal word, including every prediction, applied accelerometer,
magnetic and S correction, AW synchronization, BA evolution, Joseph gain and
attitude reset.  Stack the root as \`x_0=(x_s,x_f)\`, where \`x_s\` contains
the slow AG/BA quotient plus the physical kernel coordinate and \`x_f\`
contains the remaining nuisance root.  Stack every fresh process/sync/noise
factor once in \`s\).  Then exactly

\`y=O_s x_s+O_f x_f+A s\`,
\`x_N=T_s x_s+T_f x_f+B s\`.                                 (MW-1)

All matrices are chronological literal matrices from shipping.  No
accelerometer subword, nominal-AW mean, reset-product approximation or
independent source box is introduced.

Eliminate the nuisance root in the full joint Gaussian record.  With
\`Sigma=A A'\`, define the whitened nuisance projector

\`P_f=proj range(Sigma^-1/2 O_f)\`

and

\`J_s=O_s' Sigma^-1/2(I-P_f)Sigma^-1/2 O_s\`.                (MW-2)

Append the fictitious kernel precision \`mu=1/c\`:

\`J_(s,mu)=J_s+mu nu nu'\`.                                  (MW-3)

For the terminal state, eliminate the same nuisance root and all data
optimally using the joint minimum-action identity.  Let \`Pi_s\` be the
known-slow-root terminal covariance and \`P_(nu,s)\` the terminal covariance
with only kernel prior precision \`mu\`.  Then the exact complete-word
diameter is

\`kappa_MW(c)=lambda_max(Pi_s^-1 P_(nu,s))\`.                (MW-4)

This is not a raw observability requirement.  Directions may be weak or
instantaneously null in accelerometer rows and still have finite diameter
because the complete word combines:
- magnetic rows at every qualified service interval;
- all accelerometer epochs, including changes of nominal force;
- S pseudo-observations and their correlations with AW/LIN;
- chronological gyro-bias transport;
- AW/BA process penalties;
- actual Joseph corrections and literal reset congruences;
- terminal forgetting/process noise.

The smoother identity gives an equivalent closed-loop form.  If
\`C_k=Cov(x_s, x_k | y_<k)\`, every applied correction contributes root loss

\`Delta J_k=C_k H_k' S_k^-1 H_k C_k'\`,                     (MW-5)

transported in the same root coordinates.  Predictions contribute forgetting
through the joint terminal conditional covariance.  Thus exact row-null
directions at one epoch are harmless unless they remain invariant under the
entire corrected chronology.

The source-uniform theorem target is now directly

\`P_(nu,s)(W,c) <= K_MW(c,r) Pi_s(W)\`                       (MW-6)

for every admissible complete MOVING word whose prefixes stay in the candidate
radius.  Equivalently

\`sup_W kappa_MW(W,c,r)<infinity\`.                          (MW-7)

The physical assumptions enter only in excluding a complete-word invariant
null sequence.  A contradiction sequence with unbounded diameter would,
after normalization, have to make simultaneously:
1. every magnetic-service corrected loss vanish;
2. every accelerometer corrected loss vanish after optimal AW/BA/LIN mimic;
3. every S corrected loss vanish;
4. every AG/BA process-action penalty vanish;
5. terminal forgetting vanish.

Because process factors are retained in the joint action, an AW/BA nuisance
that changes between accelerometer epochs to maintain a force degeneracy pays
its literal process action.  Because S rows are retained, an AW/LIN mimic that
hides in accelerometer rows must also remain compatible with the integrated
chain.  This is exactly the coupling discarded by the two-epoch Schur
argument.

A quantitative proof should therefore lower-bound the **sum** of the five
nonnegative complete-word actions above, not seek a positive floor for any
single sensor family.  The already proved nuisance covariance bounds,
S-chain identity, Lemma T, magnetic service and process floors are valid
supplies for this sum.

If (MW-6) closes, no separate G0 or two-epoch \`a_2\` is needed for O1:
\`K(c,r)=K_MW(c,r)\`.  The same complete-word joint factorization gives the
known-root scalar kernel return \`D(c,r)\`; actual-gain finite-error
composition then supplies \`E(c,r)\`.  The final invariant test remains

\`D(c,r)<=c\`,
\`[1-sqrt(1-1/K_MW(c,r))]r>E(c,r)\`.

This section changes the proof target only; it does not claim a uniform
complete-word bound yet.

## Complete-word zero-action nullspace

Consider a sequence of normalized slow/root directions and nuisance/source
mimics for complete regular MOVING words whose total joint action tends to
zero.  Pass to a convergent subsequence on a fixed combinatorial word type.
All statements below are in the limiting homogeneous auxiliary system.

**1. Fresh process action rigidifies the mimic.**  Every fresh source factor
has positive covariance on its driven subspace.  Zero source action therefore
sets every AG, LIN/AW and BA fresh factor to zero.  AW sync factors also vanish.
Consequently the nuisance trajectory is no longer independently selectable at
successive observations: it is the deterministic propagation of one root
nuisance vector through the literal transitions.  In particular BA follows
its homogeneous OU decay and the complete (v,p,S,a_w) chain follows its
homogeneous linear dynamics.

**2. S loss removes the free LIN chain.**  Applied S pseudo-observations have
positive R_S.  Zero corrected S loss gives S(t_j)=0 at every applied S row.
The regular scheduler supplies the three separated S rows used in the
nuisance proof.  With zero fresh LIN/AW sources, the exact integrated-chain
Vandermonde/S-chain identity implies the only homogeneous LIN trajectory
compatible with all these zero S values and zero terminal forgetting is the
zero (v,p,S,a_w) trajectory.  Thus no AW root remains available to retune the
accelerometer force independently across epochs.

**3. Magnetic loss restricts attitude/gyro chronology.**  Zero magnetic loss
at every applied informative magnetic row gives
\`[B_k]x theta_k=0\`; hence \`theta_k=lambda_k B_k\` at those rows.  With zero
AG fresh source action, \`(theta,b_g)\` follows the literal deterministic
prediction/reset transport.  MAGNETIC SERVICE bounds service gaps.  Lemma T
then implies that any nonzero gyro-bias component that would move the
attitude away from the transported field-axis line is impossible.  Therefore
the limiting AG trajectory is the magnetic-compatible deterministic class:
gyro bias has no observable transverse/axial component left, and attitude is
the transported field-axis mode.

**4. Accelerometer loss fixes BA compatibility.**  Since the homogeneous AW
trajectory was killed in step 2, zero accelerometer loss at every applied row
reduces to

\`J_att,k theta_k + J_ba,k b_a,k=0\`.

BA has one root vector and deterministic OU decay from step 1; it cannot be
chosen separately at each epoch.  On the magnetic-compatible attitude class,
these equations define the common tilt/BA compatibility line.  The literal
shipping kernel is

\`nu=(theta_hat,0,...,0,-J_att theta_hat)\`

at a word root.  Chronological transport maps this line to the corresponding
compatibility line at subsequent rows.  Hence the complete zero-action
trajectory lies in \`span(nu)\`.

**5. Terminal forgetting excludes any hidden terminal mode.**  Zero terminal
forgetting means a remaining deterministic root direction must also survive
to the terminal state without fresh-process separation.  Steps 1--4 leave
only \`span(nu)\`; no additional nuisance terminal mode remains.

Therefore the nullspace of the complete-word joint action is contained in
the physical tilt/BA kernel:

\`Null(Action_MW) subset span(nu)\`.                           (ZN-1)

Conversely the quiet compatible construction shows why the kernel must be
retained rather than declared observable; after adding fictitious precision
\`mu nu nu'\`, the augmented action has trivial nullspace.

**Qualification still required.**  Step 2 needs four, not three, regular S
rows for the full homogeneous (v,p,S,a_w) chain.  The three-row determinant in
the nuisance proof cancels only the neutral (v,p,S) polynomial; with zero AW
process a nonzero homogeneous OU root adds an exponential mode.  Four distinct
S times give the scalar basis {1,t,t^2,psi_tau(t)}, where psi_tau is the
integrated OU contribution.  The scheduler supplies such rows, but a uniform
nonzero determinant over h/tau and the allowed S-gap intervals must be proved.
  Step 4 likewise requires that the deterministic BA
decay and transported magnetic-compatible attitude line have a common
accelerometer compatibility intersection of dimension at most one for every
admissible MOVING word.  This is a multi-epoch statement and is not implied by
the failed two-epoch force-separation lemma.  Until these two qualifications
are discharged, (ZN-1) is a proof skeleton, not a closed theorem.

## Four-S injectivity lemma for the homogeneous LIN/AW chain

For one axis with zero fresh LIN/AW sources, let the root AW mode be
\`a(t)=a_0 exp(-t/tau)\`.  Integrating through v, p and S gives

\`S(t)=S_0+p_0 t+v_0 t^2/2+a_0 psi_tau(t)\`,

where, modulo the polynomial terms already represented by \`v_0,p_0,S_0\`,

\`psi_tau(t)=tau^3[(t/tau)^2/2-t/tau+1-exp(-t/tau)]\`.       (S4-1)

Its derivatives satisfy

\`psi_tau'''(t)=exp(-t/tau)>0\`.                              (S4-2)

Hence the ordered functions
\`{1,t,t^2,psi_tau(t)}\` form a strict extended Chebyshev system on every
finite time interval for every finite \`tau>0\`.  Equivalently, for any
\`t_0<t_1<t_2<t_3\`,

\`det [1,t_i,t_i^2,psi_tau(t_i)]_(i=0..3) != 0\`,             (S4-3)

with fixed sign.  This follows directly from the generalized mean-value
formula for divided differences: after removing the Vandermonde factor, the
last divided difference equals \`psi_tau'''(xi)/3!\` for some
\`xi in (t_0,t_3)\`, and is strictly positive.

For regular A21 the S scheduler has a positive minimum separation for the
selected rows and a finite maximum gap; the shipping tuner range has
\`tau in [0.02,12]\` s.  Normalize \`t_0=0\`.  The admissible tuple
\`(t_1,t_2,t_3,tau)\` therefore lies in a compact set with
\`0<t_1<t_2<t_3\` uniformly separated.  The determinant in (S4-3) is
continuous and nowhere zero on that compact set.  Consequently

\`delta_S4 :=
 inf |det [1,t_i,t_i^2,psi_tau(t_i)]| >0\`.                  (S4-4)

Thus four selected applied S observations make the zero-source homogeneous
\`(v,p,S,a_w)\` root injective, uniformly over all regular scheduler timings
and allowed tau.  If all four S corrected losses vanish, the homogeneous
LIN/AW root is zero.  This closes qualification 1 of the complete-word
zero-action nullspace proof in real arithmetic.  A numerical value of
\`delta_S4\` is not required for qualitative coercivity; it will be needed
later for an explicit K_MW modulus.

## Multi-epoch magnetic-compatible attitude/BA intersection

Assume the zero-action conclusions already proved: all fresh AG/LIN/AW/BA
sources vanish, the homogeneous LIN/AW root is zero by the four-S lemma, and
magnetic service/Lemma T leaves no nonzero gyro-bias component.  Let
\`F_k\` be the literal invertible attitude transport from the word root to
accelerometer epoch k, including every prediction and reset.  Then

\`theta_k=F_k theta_0\`,
\`b_a,k=phi_b(t_k) R_(ba,k) b_a,0\`,                         (AB-1)

where \`R_(ba,k)\` is the literal orthogonal frame transport of the stored BA
coordinate.

Zero magnetic loss at every applied magnetic row j gives

\`F_j theta_0 in span(B_j)\`.                                (AB-2)

Thus the admissible root attitude space is the intersection

\`K_M=intersection_j F_j^-1 span(B_j)\`.                     (AB-3)

It has dimension at most one unless all pulled-back field lines coincide.  If
two do not coincide, \`K_M={0}\` and the desired result is immediate (BA is
then zero from any accelerometer row because its BA row is invertible).

In the only nontrivial case, \`K_M=span(theta_hat_0)\`.  Write
\`theta_0=lambda theta_hat_0\`.  Zero accelerometer loss at every applied
epoch k, with AW/LIN already zero, is

\`J_att,k F_k theta_0
 +R_(ba,k) phi_b(t_k)b_a,0=0\`.                              (AB-4)

Since \`R_(ba,k)\` is invertible,

\`b_a,0=-lambda phi_b(t_k)^-1
 R_(ba,k)' J_att,k F_k theta_hat_0\`                         (AB-5)

for every k.  Therefore a nonzero solution exists iff the vectors

\`q_k:=phi_b(t_k)^-1 R_(ba,k)' J_att,k F_k theta_hat_0\`     (AB-6)

are identical at all accelerometer epochs.  When they are identical, the
solution set is exactly the one-dimensional line

\`(theta_0,b_a,0)=lambda(theta_hat_0,-q)\`.                   (AB-7)

When they are not identical, only \`lambda=0\` is possible, and then
\`b_a,0=0\`.

Hence, without any force-separation premise,

\`dim{(theta_0,b_a,0): all zero magnetic and acc losses}<=1\`. (AB-8)

This is purely linear algebra: magnetic rows first reduce the root attitude
space to dimension <=1; the invertible BA row at one accelerometer epoch then
determines the entire three-component BA root from that single attitude
amplitude; additional accelerometer epochs can only remove the line, never
increase its dimension.  No assumption that the surviving line equals the
quiet physical kernel is needed for the dimension statement.

For Corollary K, however, the line used for the fictitious prior must contain
the actual complete-word nullspace.  If \`K_M=span(theta_hat_0)\` and (AB-6)
is constant, define the word's physical compatibility vector

\`nu_W=(theta_hat_0,0,...,0,-q)\`.                            (AB-9)

The complete-word nullspace is contained in \`span(nu_W)\`.  In the quiet
case this reduces to the previously stated
\`nu_ba=-J_att theta_hat\`.  On a moving word the compatible line may rotate
with the literal chronology; the recurring kernel ceiling and next-word
transport must therefore use \`nu_W\`, not assume a fixed instantaneous
quiet-kernel formula unless their equality is separately proved.

This closes the qualitative dimension-at-most-one lemma.  It also exposes a
bookkeeping correction: the complete-word kernel is the common
magnetic/accelerometer compatibility line of that word, which may be trivial;
its identification with the earlier root field-axis tilt/BA vector remains
to be checked before applying the existing scalar kernel ceiling unchanged.

## Uniform scalar treatment for the word-dependent compatibility kernel

Corollary K is basis-free: its proof uses only a rank-one root information
term \`mu nu nu'\`.  Therefore it applies unchanged to any nonzero
word-dependent null vector \`nu_W\`; no identification with the quiet kernel
is required.  Normalize the attitude component of a nontrivial compatibility
line to \`|theta_hat_W|=1\` and write

\`nu_W=(theta_hat_W,0,...,0,-q_W)\`.                          (WK-1)

From any applied accelerometer epoch k on the compatibility line,

\`q_W=phi_b(t_k)^-1 R_(ba,k)' J_att,k F_k theta_hat_W\`.     (WK-2)

The literal attitude transport is a product of orthogonal predictions and
reset factors with \`sigma_min(N)>=1\`.  For the inverse transport used in
(WK-2), \`||F_k^-1||<=1\`; forward \`||F_k||\` is finite on a candidate
retained word.  More directly, the shipping accelerometer Jacobian satisfies

\`||J_att,k||<=|f_hat,k|\`.

The estimator state clamps/tuner bounds and candidate storage give a finite
radius-local force ceiling \`F_hat,max(r)\`; BA decay obeys
\`phi_b(t)>=exp(-T_W/tau_b)\`.  Hence

\`|q_W| <= q_max(r)
 := exp(T_W/tau_b) F_hat,max(r) F_max(r)\`,                  (WK-3)

where \`F_max(r)=sup_k||F_k theta_hat_W||\` over the compact retained
coefficient class.  This is finite; it need not be small.

The existing BA marginal gives

\`P_ba,ba<=I/1600\`.

For any root covariance and normalized \`nu_W\`, block Cauchy--Schwarz yields

\`nu_W'P nu_W
 <=( sqrt(theta_hat_W'P_tt theta_hat_W)+|q_W|/40 )^2\`.       (WK-4)

Thus a uniform attitude scalar ceiling \`tau_theta(r)\` on the compact
retained root class implies

\`nu_W'P nu_W <=
 c_W(r):=(sqrt(tau_theta(r))+q_max(r)/40)^2\`.                (WK-5)

A separate pre-existing quiet-axis tilt ceiling is not required.  The
attitude ceiling may be obtained from the complete-word quotient covariance
bound itself: once the augmented complete-word action is positive on the
normalized compact kernel family, the known-root covariance plus the finite
kernel prior \`1/c\` gives a finite attitude marginal.  To avoid circularity,
use the construction/startup root ceiling for the first recurring word and
propagate the scalar ceiling with the exact kernel return \`D(c,r)\`;
Corollary K then applies word by word with each current \`nu_W\`.

Equivalently formulate the recurring invariant as

\`sup_(nu in N(r)) nu'P nu <= c\`,                            (WK-6)

where \`N(r)\` is the compact family of normalized complete-word compatibility
vectors satisfying (WK-1)--(WK-3).  The rank-one theorem is uniform over this
family whenever

\`sup_(W,nu_W in N(r)) kappa_(nu_W)(1/c)<infinity\`

and the scalar return satisfies

\`sup_W nu_(W+1)' P_(nu_W,W) nu_(W+1)<=c\`.                 (WK-7)

This is exactly the previous O1/O2 structure with a compact kernel family
instead of one fixed formula.  The BA marginal and (WK-3) make the family
bounded; continuity of the literal word maps makes it compact after including
the zero-kernel case separately.

**Result.**  Identification of \`nu_W\` with the instantaneous quiet kernel is
unnecessary.  The existing rank-one Corollary K generalizes without change,
and the scalar invariant generalizes to the compact family \`N(r)\`.  What
remains quantitative is the same complete-word diameter/return bound
\`K_MW(c,r),D(c,r)\`; there is no new observability obstruction caused by the
word-dependent kernel direction.

## Pointwise nullspace does not imply uniform complete-word coercivity

The earlier claim that trivial augmented nullspace plus compactness implies a
uniform `g_MW(c,r)>0` is RETRACTED.  The compatibility/kernel line can change
rank or become tangent along an admissible sequence.  Pointwise injectivity
modulo the word kernel does not imply a uniform positive nonzero singular
value across such limits.

The later finite-horizon detectability formulation is controlling.  What is
needed is a uniform relative inequality between terminal persistence and
observation action, not a uniform Euclidean information floor.

## Four-S injectivity: quantitative bound retracted

The qualitative four-S injectivity lemma remains valid because
`psi_tau'''(t)=exp(-t/tau)>0` for every finite `tau>0`, so the generalized
Vandermonde determinant is nonzero at four distinct S times.

The previously claimed explicit determinant floor is RETRACTED.  The proof
incorrectly used `tau<=12` to infer
`exp(-t/tau)>=exp(-t/12)`; the inequality is reversed.  Since the admitted
range includes `tau_min=0.02`, the divided-difference lower bound based only
on `min psi'''` is exponentially tiny and the displayed >6.82e4 claim is
false.  For example at times (0,8,16,24) and `tau=.02`, the determinant is
approximately `.008192(1-exp(-400))^3`, about .008192.

Thus four-S zero-source injectivity is qualitative only at this stage.  No
explicit four-S singular-value modulus is currently certified.

## Moving-kernel quotient Riccati diameter

For each complete word W retain its normalized compatibility kernel nu_W.
Do not require a uniform positive nonzero eigenvalue of the data information.
With kernel precision mu=1/c, let Pi_W be the known-root terminal covariance
and P_nu,W the terminal covariance with only the rank-one kernel prior. Define

`kappa_Q(W,c)=lambda_max(Pi_W^-1 P_nu,W)`.

This is basis-free. In root coordinates split the moving kernel line from its
orthogonal quotient and write

`J=[[j_nn,j_nq'];[j_nq,J_qq]]`.

After adding mu on the kernel line, block elimination gives the quotient Schur
information

`S_q=J_qq-j_nq j_nq'/(j_nn+mu)`.

Write the terminal excess map in Pi-whitened coordinates as
`H=[h_n,H_q]`. The corresponding quotient terminal map is

`H_eff=H_q-h_n j_nq'/(j_nn+mu)`.

The dangerous quantity is the relative pair (H_eff,S_q), not
`lambda_min^+(J)`. A quotient information eigenvalue may tend to zero
without making the Riccati diameter diverge when the terminal excess vanishes
at the same rate.

The exact target is the Loewner relative-action inequality

`H_eff' H_eff <= (K_Q-1) S_q`,

together with the scalar kernel block bound
`|h_n|^2/(j_nn+mu)<=K_n-1`. Equivalently, for the optimal complete-word
minimum-action reader L,

`(T-LA)(T-LA)' <= (K_Q-1) Pi_W`.

These formulations remain meaningful at rank-changing words by shorted
operator/Moore-Penrose continuation.

A finite continuous extension through a quotient rank loss requires

`Null(S_q) subset Null(H_eff)`.                             (RQ-1)

This is a range condition, not a principal-angle condition. It says that a
root quotient direction carrying zero complete-word reduced information also
carries zero terminal excess relative to the known-root covariance. The
zero-action theorem supplies the data-null classification; what remains is to
show that the same zero-action nuisance/source mimic reproduces the terminal
image, establishing (RQ-1) for the literal joint operator.

If (RQ-1) holds at every limiting word, the generalized relative-action ratio
extends continuously over the compact moving-kernel word family. Its maximum
is then finite, giving `K_MW(c,r)<infinity` without any uniform
`sigma_min^+` or principal-angle floor.

## Exact range/shorting test for the moving quotient

The joint minimum-action identity decides (RQ-1) exactly.  After nuisance and
kernel-line elimination, write the reduced auxiliary model as

`y=O_q q+A s`,
`x_N=T_q q+T s`,
`Sigma=A A'>0`.

Define

`J_q=O_q' Sigma^-1 O_q`,
`Ttilde_q=T_q-T A' Sigma^-1 O_q`.                          (RS-1)

The known-root terminal covariance is

`Pi=T(I-A'Sigma^-1 A)T'`,

and the diffuse quotient terminal covariance is

`Pi+Ttilde_q J_q^dagger Ttilde_q'`

whenever the Moore-Penrose expression is finite on the information range.
For any v in Null(J_q), positivity of Sigma gives

`O_q v=0`.                                                  (RS-2)

The desired range condition is therefore exactly

`O_q v=0 => Ttilde_q v=0`.                                 (RS-3)

But (RS-3) is **not an algebraic consequence** of the joint Gaussian model.
If `O_q v=0`, then the source least-squares mimic is zero and
`Ttilde_q v=T_q v`.  Thus (RS-3) reduces to

`Null(O_q) subset Null(T_q)`.                              (RS-4)

In words: every complete-word data-null quotient root direction must also
have zero deterministic terminal image.  A root mode can be invisible to all
measurements yet survive to the terminal state; process noise in Pi does not
make its deterministic image disappear.  This is precisely why information-
only contraction and raw observability were insufficient earlier.

For the literal complete-word zero-action classification, data-null directions
are contained in the one-dimensional compatibility family before the kernel
prior.  After quotienting by the **exact same** compatibility line nu_W, the
qualitative result gives Null(O_q)={0} for each fixed nondegenerate word.
At a rank-changing limit, however, a new quotient data-null direction may
appear as the compatibility line changes/disappears.  To obtain a continuous
diameter extension one must show its terminal image tends to zero at the same
rate; pointwise nullity does not establish that.

Therefore (RQ-1) is not yet closed.  It is equivalent to a quantitative
forgetting statement for the near-null quotient direction, not merely to the
zero-action nullspace theorem.  The existing complete-word machinery contains
the right quantity: the pair `(O_q,T_q)` or, equivalently, the generalized
ratio

`sup_v |Ttilde_q v|_Pi^-1^2 / |O_q v|_Sigma^-1^2`.         (RS-5)

A finite bound on (RS-5) is exactly the kernel-bounded Riccati diameter sought;
using (RQ-1) to prove it would be circular unless (RS-4) is established
uniformly at limiting words by an independent terminal-retention argument.

This closes the proposed automatic shorting route as a dead end.  The next
valid target is to bound the generalized pair (Ttilde_q,O_q) directly using
the literal prediction/process structure -- a detectability/finite-horizon
observability inequality -- rather than seeking a separate information floor
or assuming the null inclusion.

## Finite-horizon detectability as regularized backward readout

For the moving quotient, exact root cancellation is unnecessarily strong near
a rank-changing word.  Use the kernel-regularized joint reader.

After nuisance elimination and kernel splitting, let the quotient auxiliary
model be

`y=O_q v+A s`,
`x_N=T_q v+T s`,
`Sigma=A A'`.

Let `J_q=O_q' Sigma^-1 O_q` and let the kernel prior induce the positive
regularization inherited from the full block Schur complement; denote it
`M_q>=0`.  The minimum-action regularized reader is

`L_q=T A' Sigma^-1
 +(T_q-T A' Sigma^-1 O_q)(J_q+M_q)^dagger O_q' Sigma^-1`.   (DET-1)

Its terminal root residual is

`R_qroot=Ttilde_q[I-(J_q+M_q)^dagger J_q]`,                (DET-2)

and its fresh-source action is the corresponding backward-readout sum.  The
rank-one kernel prior makes the full reader finite even as the moving quotient
changes rank.

The desired detectability inequality

`|Ttilde_q v|_Pi^-1^2 <= C_det |O_q v|_Sigma^-1^2`         (DET-3)

is equivalent, on Range(J_q), to

`Ttilde_q' Pi^-1 Ttilde_q <= C_det J_q`.                   (DET-4)

For a fixed word its sharp constant is

`C_det(W)=lambda_max(J_q^dagger/2
 Ttilde_q'Pi^-1 Ttilde_q J_q^dagger/2)`,                   (DET-5)

with infinity if Null(J_q) is not contained in Null(Ttilde_q).  Thus
`K_MW<=1+C_det`.

Backward readout supplies a constructive upper bound for (DET-5).  Initialize
a terminal quotient residual Y at the end of the word and run backward:
prediction adds `(YU)(YU)'` and maps `Y<-YF`; observation block L_i adds
`(L_i V_i)(L_i V_i)'` and maps `Y<-Y-L_iH_i`; reset maps
`Y<-YG_i`.  Instead of requiring exact root cancellation, choose the
blocks L_i by the regularized normal equations associated with
`J_q+M_q`.  The resulting action is exactly the numerator represented in
(DET-5), with all shared process factors and literal resets retained.

This formulation shows what must be proved for a source-uniform constant:
the regularized backward-reader action divided by the root observation action
must remain bounded as a quotient direction approaches the moving kernel.
No absolute observation singular-value floor is required.

However the existing zero-action/dimension theorem is insufficient to bound
this ratio.  It gives pointwise injectivity modulo the kernel but no rate at
which terminal persistence vanishes relative to observation action near a
rank-changing word.  The regularized reader prevents an algebraic blow-up in
the construction, but (DET-3) itself can still fail if a near-null direction
has O(epsilon) observation amplitude and O(1) terminal image.

Therefore finite-horizon detectability is now the exact remaining quantitative
property; it has not been proved by the prior compactness arguments.  A valid
next step must exploit the literal terminal dynamics to show that any
compatibility-line rotation producing O(epsilon) observation action also
produces O(epsilon) terminal quotient image.  Without such a terminal
retention estimate, assigning a finite C_det would be circular.

## Terminal retention from compatibility-line variation

Let K_k be the one-dimensional compatibility line pulled back to the root from
the corrected history through epoch k, and let p_k be its orthogonal projector
in the normalized slow root metric.  Exact prediction/reset transport carries
a vector aligned with K_k without creating a quotient component; quotient
persistence is created only when the compatible line changes.

For a normalized root vector v decompose recursively

`v=p_k v+(I-p_k)v`.

At the next constraint block,

`(I-p_(k+1))p_k v=(p_k-p_(k+1))p_k v`.                    (TR-1)

Thus the terminal quotient image obeys the exact telescoping form

`T_q v=sum_k T_(N<-k) E_k (p_k-p_(k+1)) p_k v
        +T_(N<-0)(I-p_0)v`,                                 (TR-2)

where E_k is the literal embedding at the kth constraint block and
T_(N<-k) is the subsequent deterministic terminal transport.  There is no
charge for motion along the current compatibility line.

Each projector change is controlled by the residual constraint operator that
defines the new line.  If C_k is the whitened magnetic/accelerometer/S block
after nuisance shorting and K_k=Null(C_<=k) has dimension one, the
Davis--Kahan/sine relation gives

`||(p_k-p_(k+1))p_k v||
 <= ||C_(k+1) p_k v|| / gap_(k+1)`,                         (TR-3)

where gap_(k+1) is the smallest positive singular value of C_(k+1) restricted
to K_k^perp.  Multiplying by terminal transport gives

`||T_q v||_(Pi^-1)
 <= A_0 ||(I-p_0)v||
   +sum_k A_k/gap_(k+1) ||C_(k+1)p_k v||`,                  (TR-4)

with
`A_k=||Pi^-1/2 T_(N<-k)E_k||`.

Cauchy--Schwarz then yields the detectability constant

`C_ret <= A_0^2/gap_0^2
 +sum_k A_k^2/gap_(k+1)^2`,                                (TR-5)

provided the block residual energies are the same terms appearing in
`|O_qv|_Sigma^-1^2`.

This is the desired terminal-retention mechanism: only compatibility-line
rotation is charged, and aligned persistence is quotiented out.

But (TR-3) exposes the same quantitative issue as the retired global
principal-angle proof, now locally: a finite C_ret requires a lower positive
restricted gap for each line-changing constraint block.  Current assumptions
allow a new constraint to become arbitrarily tangent to the existing
compatibility line.  In that case both the line rotation and residual are
small, but their ratio is governed by the local restricted singular value,
which has no proved lower bound.

There is no generic improvement from telescoping: for a two-dimensional
example C_epsilon=[epsilon,0], the compatible line is fixed while observation
energy is epsilon^2 and a later deterministic terminal map can retain the
first coordinate with O(1) amplitude.  Exact transport and compactness do not
bound the ratio.

Therefore a source-uniform C_ret is not established by compatibility-line
variation under the current assumptions.  The terminal-retention route is
equivalent to a local transversality/detectability modulus.  This is not a
new algebraic gap: it is the same near-null/terminal-persistent obstruction
identified by (DET-3).

The proof must now decide between two genuinely different mechanisms:
(1) derive a restricted gap from an already existing shipping invariant or
scheduler/service condition not yet used, or
(2) accept that the current physical assumptions permit arbitrarily weak
detectability and cannot yield a source-uniform contraction margin.  Numerical
carried-word gaps cannot decide this theorem question.

## Continuous soft-kernel Riccati diameter

The exact-kernel quotient is discontinuous when the compatibility nullspace
changes dimension. Avoid that quotient. Carry instead a continuous unit
direction n(W) selected from the least-information tilt/BA eigenspace, with
sign fixed by the transported magnetic field, and append finite precision
`mu n n'` whether or not the data nullspace is exactly one-dimensional.

Define

`J_mu(W)=J(W)+mu n(W)n(W)'`,
`P_mu(W)=Pi+Phi_tilde J_mu(W)^-1 Phi_tilde'`,
`kappa_soft(W)=lambda_max(Pi^-1 P_mu(W))`.

No direction is abruptly promoted from kernel to quotient when an exact
kernel disappears. If a formerly exact kernel acquires information epsilon^2,
the same finite prior mu remains on that direction and its contribution is
bounded by `1/(mu+epsilon^2)` rather than `1/epsilon^2`.

Let E span n^perp and block
`J=[[j_nn,j_nq'];[j_nq,J_qq]]`. The exact Schur complement is

`S_soft=J_qq-j_nq j_nq'/(j_nn+mu)`.

Because `j_nn+mu>=mu`, rank changes of J along n are harmless. The remaining
requirement is coercivity only on directions uniformly separated from the
chosen soft line. This is weaker than the retired exact-kernel principal-angle
floor: when a second weak direction approaches n, the least-information
eigenvector n rotates continuously with the weak subspace rather than
disappearing.

A basis-free formulation uses the two smallest eigenvalues of the tilt/BA
information block. If lambda_2(W), the second eigenvalue after the soft
direction, has a source-uniform positive floor, then

`S_soft >= lambda_2,bar I`

on the soft quotient and

`kappa_soft <= 1 + ||Pi^-1/2 Phi_tilde||^2/lambda_2,bar`

(with the sharper block reader retained for actual constants).

The qualitative complete-word nullspace theorem already proves nullity at
most one for each word. What it did not prove was a uniform smallest positive
eigenvalue across rank-changing words. For the soft-kernel theorem the
relevant compactness quantity is instead lambda_2, counting eigenvalues from
zero with multiplicity. Unlike `lambda_min^+`, lambda_2 is continuous through
a one-dimensional kernel appearing/disappearing. If every limiting word has
nullity at most one, compactness now legitimately implies

`inf_W lambda_2(W)>0`,

provided the admissible retained word class itself is compact and the
complete-word information matrix J(W) is continuous on that closed class.

This repairs the earlier compactness error: continuity applies to the fixed
ordered eigenvalue lambda_2, not to `lambda_min^+`. It also eliminates the
kernel-disappearance epsilon^-2 counterexample because finite precision is
kept on the continuously selected soft direction on both sides of the rank
change.

The next obligations are therefore (i) verify compactness/closedness of the
radius-local complete-word coefficient class including scheduler/event-type
limits; (ii) verify continuity of J on each event stratum and across allowed
event-boundary limits, or use finitely many closed strata; and (iii) define a
continuous/measurable soft direction n(W). A globally continuous eigenvector
is unnecessary for the bound: the rank-one projector onto the least
eigenspace suffices when lambda_1<lambda_2; at equality any minimizing
projector gives the same conservative two-dimensional treatment.

If these topological obligations close, a finite source-uniform soft-kernel
diameter follows without AW-gain entry, determinant sign crossing, or
finite-horizon detectability of the disappearing exact kernel.
## Finite closed event strata for regular A21 words

Fix the regular A21 word horizon T_W and the qualified sample interval
`h in [h_min,h_max]` with h_min>0. Then the number of prediction slots is
bounded by `N_max=ceil(T_W/h_min)+1`. At each slot the shipping operation
alphabet is finite: prediction; accelerometer attempted/applied or rejected;
S due/applied or failed-safe; magnetic attempted/applied or rejected; BA
mode/projection branch; AW sync due/not due; and the finite reset/projection
branches. Hence there are finitely many combinatorial event words on T_W.

For a fixed combinatorial word sigma, collect continuous coefficients
(dt values, states, covariance, physical/tuner parameters, measurements,
clocks and references) into z. Literal prediction, Joseph correction, reset
and covariance-sync maps are continuous wherever their selected branch has
strictly positive LDLT pivots and finite inputs. Define a closed stratum by
replacing each branch predicate with a closed assignment convention at
equality (for example the code's actual <= versus < choice). Gate-boundary
points that select the opposite code branch are placed in that neighboring
stratum rather than duplicated.

The closure of one fixed applied-event formula need not equal a shipping
stratum: as an LDLT pivot approaches the fail-safe boundary the applied
formula may cease to be the code path. Therefore continuity is asserted only
on each code-selected closed branch, not across a branch switch. The finite
union of these branch strata covers all regular words.

MAGNETIC SERVICE is evaluated from actually applied informative events. For
a fixed applied-magnetic pattern its service Gram is continuous in z and the
condition `lambda_min(sum G_k'G_k)>=mu_M` is closed. A limit in which an
event becomes rejected belongs to the rejected neighboring stratum; that
stratum is admissible only if its remaining actually applied events still
satisfy the same closed service inequality. Thus service prevents loss of the
last required magnetic information at a boundary unless other applied events
already retain the floor.

Do not stratify by exact nuisance rank: exact-rank sets are generally not
closed. Instead retain the raw fixed-dimensional joint factor action
`Q_W(x)=min_z ||A_W x+D_W z||^2`, with nuisance/process action rows included
in `D_W`. Structural zero-action directions are moved into the slow root once
and for all. The proved action floors make the penalized nuisance operator
coercive, so minimizers are bounded and `Q_W` is continuous/lower-
semicontinuous through observation-rank changes. The symmetric slow
information `J(W)` is defined by this quadratic form, without a data-dependent
pseudoinverse.

The remaining topological requirement is compactness of each substratum.
Bounds on dt, tuner clamps, bias/state retained radius, covariance ceilings,
physical/reference bounds and fixed horizon make the continuous coefficient
set bounded. Closed branch predicates, closed MARINE MOTION/IMU BIAS bounds
and closed MAGNETIC SERVICE make it closed, provided the all-time physical
continuation variables are represented by the already assumed compact
finite-word trace class. That last trace compactness is not automatic from
bounded acceleration/jerk alone and must use Arzela--Ascoli plus the bounded
potential/velocity/displacement conditions.

If compactness holds and the complete-word nullity<=1 theorem applies on
every boundary substratum, ordered eigenvalue continuity gives
`lambda_2>0` pointwise and therefore a positive minimum on each substratum;
the finite minimum over strata is positive.

Boundary nullity is the remaining substantive issue. The earlier zero-action
proof uses four applied S observations and recurring magnetic service. A
boundary stratum may lose an S application through an LDLT fail-safe; unlike
MAGNETIC SERVICE, there is currently no first-class S-service assumption
requiring four actually applied S rows. Therefore the existing nullity<=1
proof does not automatically extend to every closure stratum. A boundary
with too few applied S corrections can retain additional homogeneous LIN/AW
null directions while still satisfying MAGNETIC SERVICE.

Scheduled S-row survival is in fact guaranteed in real arithmetic on the
retained shipping domain. The S innovation covariance is
`P_SS+R_S`; covariance PSD gives `P_SS>=0`, while the shipping tuner clamp
`sigma_S>=1e-6` gives `R_S>=1e-12 I`. Hence `P_SS+R_S` is uniformly SPD and
the mathematical LDLT cannot fail. Therefore every scheduled regular S update
is applied, and the four-S zero-action argument extends to event-boundary
strata without a new S-service assumption.

The remaining boundary question is then magnetic/accelerometer compatibility.
MAGNETIC SERVICE is closed and stated using actually applied magnetic rows,
so every admissible boundary stratum retains its magnetic information floor.
Accelerometer invalid-input rejection is excluded on the finite physical
retained class by finite measurements/temperature; its LDLT likewise has
`S_acc>=Racc>=1e-8 I` from the sigma_acc clamp in real arithmetic. Thus
scheduled finite accelerometer updates are also mathematically applied.

Accordingly the operation pattern relevant to the nullity theorem is fixed by
the regular scheduler except for magnetic acceptance/service patterns, of
which there are finitely many and each admissible closed pattern retains
MAGNETIC SERVICE. Under the existing complete-word zero-action proof this
supports nullity<=1 on every admissible stratum closure.

## Compactness of finite-word physical traces

Fix a complete MOVING A21 word interval I=[0,T_W]. The existing MARINE
MOTION bounds give uniform pointwise bounds on p,v,a and the a.e. jerk bound
`|a_dot|<=J_max`. Thus a is equibounded and J_max-Lipschitz. Arzela--Ascoli
gives a uniformly convergent subsequence of a. Bounded initial v and p plus
`v(t)=v(0)+int a` and `p(t)=p(0)+int v` then give uniform convergence of v
and p. The kinematic identities pass to the limit.

The jerk condition is closed because a uniform limit of J_max-Lipschitz
functions is J_max-Lipschitz. Pointwise p/v/a bounds are closed. The
bounded-potential condition `|int_t1^t2 p|<=P_AC` is closed under uniform
convergence of p for every pair t1<t2.

Physical accelerometer/gyro bias traces are uniformly bounded with uniform
derivative bounds, hence are equibounded/equi-Lipschitz and compact in C0 by
Arzela--Ascoli; their derivative bounds remain valid in the limit.

For attitude use R(t) in SO(3). SO(3) is compact, and the nominal/physical
angular-rate bound makes R equi-Lipschitz through `R_dot=R[omega]x` in the
integral sense. A uniformly convergent subsequence therefore has a limit in
SO(3) satisfying the same attitude kinematics. The additional bounded rate
traces needed by lever-arm coefficients are treated by their existing
regularity bounds.

The MOVING excitation condition is closed: on each complete T_E subwindow,
attitude span is the maximum of a continuous rotation-distance function over
a compact pair-time set. Uniform convergence of R preserves
`span>=theta_E`. STILL/transition histories belong to their separate regime
strata and are not inserted into this recurring MOVING word class.

All-time extendability is also closed. Take a sequence of globally admissible
continuations whose restrictions converge on I. Apply the same uniform
bounds and Arzela--Ascoli on [-n,n] for n=1,2,... and use a diagonal
subsequence. The local limits agree on overlaps and define a global
continuation. The derivative, pointwise and bounded-potential inequalities
pass to that limit. Hence restrictions to I of globally admissible histories
form a closed compact trace class.

Sensor/frontend samples on the finite scheduler are continuous evaluations
of these traces and bounded finite-dimensional frontend/reference variables.
Combining this physical compactness with the already bounded closed
estimator/covariance/tuner/clock variables proves compactness of every closed
event/rank stratum described above.

Therefore, assuming the established boundary nullity<=1 result on every
admissible stratum, the ordered second eigenvalue lambda_2(J(W)) is continuous
and strictly positive on each compact stratum. It attains a positive minimum
there. There are finitely many strata, so

`lambda2_bar := min_strata min_W lambda_2(J(W)) > 0`.

This is a legitimate source-uniform **existence** modulus for the soft-kernel
route. No numerical value is asserted. The argument uses the fixed ordered
eigenvalue lambda_2, which is continuous through a one-dimensional rank
change, rather than the discontinuous `lambda_min^+` used in the retracted
coercivity argument.
## Soft-kernel Riccati diameter from the ordered spectral gap

For each complete word W let the ordered eigenvalues of the nuisance-reduced
slow information be `0<=lambda_1<=lambda_2<=...`, and choose any unit
eigenvector n_W for lambda_1. Put `mu=1/c` and

`J_soft=J+mu n_W n_W'`.

Because n_W is an eigenvector of J, the eigenvalues of J_soft are exactly
`lambda_1+mu, lambda_2,...`. Therefore

`J_soft >= m(c,r) I`,
`m(c,r):=min(mu,lambda2_bar)=min(1/c,lambda2_bar)>0`.        (SK-1)

No continuity of the chosen eigenvector is needed for this inequality. At a
multiple least eigenvalue, any minimizing eigenvector gives the same lower
bound; lambda2_bar>0 prevents a two-dimensional zero eigenspace.

The exact Riccati identity gives

`P_soft=Pi+Phi_tilde J_soft^-1 Phi_tilde'`.

Hence

`Pi^-1/2 P_soft Pi^-1/2
 <= I + (1/m) Pi^-1/2 Phi_tilde Phi_tilde' Pi^-1/2`.

Let
`Hbar(r):=sup_W ||Pi_W^-1/2 Phi_tilde_W||_2^2`
over the compact retained word class. Pi has the already proved regular
known-root positive process floor and all literal word maps are continuous,
so Hbar(r)<infinity. Consequently

`K_soft(c,r):=1+Hbar(r)/min(1/c,lambda2_bar)`               (SK-2)

is a valid source-uniform soft-kernel Riccati diameter:

`P_soft(W,c)<=K_soft(c,r) Pi_W`.

This is an existence-level explicit formula in two compactness constants
`lambda2_bar` and `Hbar`; neither has yet been numerically extracted.

### Soft scalar return

The invariant scalar must use the same least-information direction family,
not the old exact physical kernel. Let n_W be the root soft direction for
word W and n_next any least-information unit direction at the next root.
Define

`d_soft(W):=n_next' Pi_W n_next`

and

`dbar_soft(r):=sup_W,n_next d_soft(W)<infinity`,             (SK-3)

where finiteness follows from compactness and continuity of Pi; no full-state
recurring covariance ceiling is required.

From `P_end<=P_soft<=K_soft Pi`,

`n_next' P_end n_next <= K_soft(c,r) dbar_soft(r)`.

Therefore the soft scalar invariant closes whenever

`D_soft(c,r):=K_soft(c,r)dbar_soft(r) <= c`.                (SK-4)

This replaces the earlier exact-kernel O2 fixed point. It remains linked to
the same word through Pi before taking the source-uniform sup.

There is an important asymptotic consequence. Since
`m=min(1/c,lambda2_bar)`, for `c>=1/lambda2_bar`,

`K_soft=1+Hbar c`.

Then `D_soft<=c` requires
`dbar_soft(1+Hbar c)<=c`, i.e.
`dbar_soft Hbar<1` and
`c>=dbar_soft/(1-dbar_soft Hbar)`.

Thus existence of finite lambda2_bar alone is not enough for the scalar
invariant: the dimensionless linked product `dbar_soft Hbar` must be <1 in
the large-c branch (or the small-c branch must close). This is now the
controlling quantitative/existence test for O2.
## Linked soft-kernel return/diameter action

Do not use the product of separate suprema dbar_soft Hbar. For one word W
and its next-root soft direction n_+, define

`d_W=n_+' Pi_W n_+`,
`H_W=Pi_W^-1/2 Phi_tilde_W`.

The large-c scalar return obtained from the coarse spectral diameter contains
`d_W ||H_W||^2`. Preserve the same word and define

`Gamma_soft(r):=sup_W d_W ||H_W||_2^2`.                    (LK-1)

This is the correct linked replacement for `dbar_soft Hbar`; always
`Gamma_soft<=dbar_soft Hbar`, often strictly.

Let `z_W=Pi_W^(1/2)n_+`, so `|z_W|^2=d_W`. Then

`d_W ||Pi_W^-1/2 Phi_tilde_W||^2
 = ||z_W||^2 ||Pi_W^-1/2 Phi_tilde_W||^2`.                 (LK-2)

A direct rank-one terminal bound is sharper for the scalar return:

`n_+' Phi_tilde J_soft^-1 Phi_tilde' n_+
 = ||J_soft^-1/2 Phi_tilde' n_+||^2`.                       (LK-3)

Using `J_soft>=m I` gives

`<= (1/m) ||Phi_tilde' n_+||^2`.                            (LK-4)

Thus the scalar invariant does not actually require the full operator H_W.
Define the same-word terminal coupling

`g_W(c):=n_+' Phi_tilde_W J_soft(W,c)^-1 Phi_tilde_W' n_+`,
`gbar_soft(c,r):=sup_W g_W(c)<infinity`.                    (LK-5)

Then exactly

`n_+' P_soft n_+ = d_W + g_W(c)`.                          (LK-6)

and the source-uniform scalar return can be taken as

`D_soft(c,r)=sup_W [d_W+g_W(c)]`.                          (LK-7)

This dominates neither factor separately and preserves the linked word.

The previously stated large-c criterion based on
`||Phi_tilde' n_+||^2<1` is RETRACTED: that Euclidean quantity is not
coordinate/units invariant and the factor `1/m` was already an upper bound,
not an exact scalar recurrence coefficient.

The invariant condition is simply

`sup_W [d_W+g_W(c)] <= c`,                                  (LK-8 corrected)

with
`g_W(c)=n_+' Phi_tilde (J+mu n_W n_W')^-1 Phi_tilde' n_+`
and `mu=1/c`. This is dimensionally and geometrically correct.

Equivalently normalize terminal directions by Pi: put
`z_+=Pi^(1/2)n_+/sqrt(d_W)` when d_W>0. Then

`g_W(c)/d_W = z_+' Pi^-1/2 Phi_tilde J_soft^-1
               Phi_tilde' Pi^-1/2 z_+`,                    (LK-9)

which is bounded by `kappa_soft-1` and is dimensionless. But requiring this
ratio `<1` is sufficient, not necessary; the exact fixed-point inequality
above is controlling.

Thus O2 should now be attacked by the exact linked scalar function
`F_W(c)=d_W+g_W(c)-c`. For each word, `g_W(c)` is continuous and
nondecreasing in c (because mu=1/c decreases and inverse information grows).
A source-uniform invariant exists if one can find c with
`sup_W F_W(c)<=0`. Compactness makes the supremum attained and continuous
under the soft-kernel stratum construction. The next task is to study its
large-c slope using the exact least-information spectral component, not a
coordinate-dependent Euclidean norm.
## Exact large-c soft-kernel persistence coefficient

Diagonalize one word's nuisance-reduced slow information as
`J=sum_i lambda_i e_i e_i'`, choose the soft direction `n_W=e_1`, and put
`mu=1/c`. Then exactly

`(J+mu n_W n_W')^-1
 =(lambda_1+mu)^-1 n_W n_W'
  +sum_(i>=2) lambda_i^-1 e_i e_i'`.                        (AL-1)

For next-root soft functional n_+, set
`a_i=n_+' Phi_tilde e_i`. Therefore

`g_W(c)=a_1^2/(lambda_1+1/c)
       +sum_(i>=2) a_i^2/lambda_i`.                         (AL-2)

If lambda_1>0, g_W(c) remains bounded and its large-c linear coefficient is
zero. If lambda_1=0, the formula is exact:

`g_W(c)=c |n_+' Phi_tilde n_W|^2 + beta_W`,
`beta_W=sum_(i>=2) a_i^2/lambda_i`.                          (AL-3)

Hence

`alpha_W=0` when lambda_1>0,
`alpha_W=|n_+' Phi_tilde n_W|^2` when lambda_1=0.            (AL-4)

There is no asymptotic remainder on an exact-kernel word.

However AL-4 must be interpreted with the same kernel-functional
normalization used to define the scalar ceiling c. Euclidean unit
eigenvectors mix attitude and BA coordinates with different physical units;
under a rescaling of root coordinates the numerical comparison
`alpha_W<1` changes. The earlier physical kernel formulation implicitly fixed
a functional nu and bounded `nu'Pnu`; the soft spectral construction must
likewise declare a root metric/dual functional normalization.

After fixing a positive root metric M, normalize the soft functional in the
corresponding dual metric and insert the matching M factors in AL-4. The
source-uniform large-c slope is then the squared induced transfer of the
current exact data-null functional into the next normalized soft functional.

The existing theorem set does not prove this normalized transfer is strictly
below one. MAGNETIC SERVICE and four-S control directions transverse to the
exact compatibility kernel; on an exact data-null line they supply no
measurement loss. BA decay helps its BA component, but the attitude component
is transported nearly neutrally. Therefore lambda2_bar>0 does not imply
`sup alpha_W<1`.

Conclusion: alpha_W is now derived exactly, but O2 remains an independent
kernel-persistence obligation. The next valid step is to restore the
physically declared kernel functional/metric (rather than arbitrary Euclidean
soft eigenvector normalization) and analyze its exact one-word transfer using
BA decay plus attitude/reset chronology.
## Physical normalization and exact compatibility-line transfer

Use the physical attitude amplitude as the scalar kernel coordinate. At a
word root with a nontrivial complete-word compatibility line choose its
generator

`nu_W=(theta_hat_W,0,...,0,-q_W)`, with `|theta_hat_W|=1`,

where theta_hat_W is the root attitude-error direction about the pulled-back
world magnetic-field axis and q_W is fixed by all accelerometer compatibility
equations. A kernel state is `x=lambda nu_W`; lambda is therefore an angle
amplitude in radians. The scalar covariance ceiling is the variance of this
declared lambda functional, not the Euclidean norm of the mixed attitude/BA
vector.

On an exact zero-loss word all fresh root-source factors are zero for the
homogeneous transfer. Four-S kills LIN/AW, magnetic service/Lemma T kills the
gyro-bias quotient, and the surviving attitude obeys the deterministic
attitude/reset transport

`theta_N=F_W theta_0`.

By definition of the next compatibility generator, if the exact line
persists across the word then

`F_W theta_hat_W = s_W theta_hat_(W+1)`                     (PT-1)

for a scalar signed attitude-amplitude transfer s_W. Because both generators
are normalized to unit attitude norm,

`|s_W|=|F_W theta_hat_W|`.                                  (PT-2)

The BA component propagates homogeneously as
`b_a,N=Phi_b,W (-lambda q_W)` (including its literal frame transport).
Exact compatibility at the next root requires

`Phi_b,W q_W = s_W q_(W+1)`.                               (PT-3)

Thus BA decay does not independently multiply the scalar lambda. Instead it
constrains which exact compatibility lines can persist: q_(W+1) must absorb
the relative BA decay compared with attitude transport.

The physically normalized large-c persistence coefficient is therefore

`alpha_W^phys = s_W^2`                                      (PT-4)

for an exact data-null word whose compatibility line maps into the next
compatibility line; it is zero in the large-c sense on words with no exact
kernel.

Now inspect the literal attitude transport. Prediction of the homogeneous
attitude error with zero gyro-bias is a rotation/orthogonal transport in the
world-frame representation. Quaternion reset uses the exact left-error
coordinate change; along a zero-innovation/zero-loss compatibility trajectory
the nominal correction associated with the homogeneous perturbation is zero,
so the homogeneous physical attitude perturbation is transported without
dissipation. In the world-frame historical coordinates used by the proof,
its norm is preserved. Hence for an exact persistent compatibility line

`|F_W theta_hat_W|=1`,

and therefore

`alpha_W^phys=1`.                                           (PT-5)

BA OU decay does not make this scalar mode contract: equation (PT-3) changes
the BA-to-tilt ratio q of the next compatibility generator, while the
physically normalized tilt amplitude remains unit persistent.

Consequently, if the admissible recurring class contains an exact
compatibility line that persists from one complete word root to the next,
then

`sup_W alpha_W^phys >=1`.

The soft-kernel O2 large-c strict-contraction condition cannot then be proved
from BA decay alone. A strict `<1` requires an additional mechanism that
attenuates the physical attitude amplitude along the exact compatibility
line, or a theorem showing no exact compatibility line can persist across
successive recurring words under the current assumptions.

This is a conditional obstruction, not yet an existence proof of a persistent
exact-compatible recurring execution. It clarifies the role of BA decay:
decay rotates/reweights the compatibility line; it does not itself dissipate
the physically normalized tilt scalar.
## Persistence of an exact tilt/BA compatibility line across words

Assume a nontrivial exact compatibility line persists through consecutive
regular A21 words and normalize every generator by unit physical attitude
amplitude. The previous transfer identity gives

`Phi_b,W q_W = s_W q_(W+1)`, with `|s_W|=1`.

Active A21 BA prediction has `0<Phi_b,W<1` over every positive-duration
complete word. Hence

`|q_(W+n)| = (prod_j Phi_b,j) |q_W|`

and, using the positive minimum word duration and finite tau_b ceiling,

`|q_(W+n)| <= phi_bar^n |q_W| ->0` for some phi_bar<1.       (PER-1)

Therefore indefinite exact persistence forces the compatibility lines toward
the zero-BA field-axis line.

At each accelerometer epoch k in word W, exact compatibility is

`R_ba,k phi_b,k q_W = J_att,k F_k theta_hat_W`              (PER-2)

up to the fixed sign convention. Thus PER-1 implies that along an indefinitely
persistent sequence the transverse nominal specific-force action on the
transported magnetic-axis tilt tends uniformly to zero at the applied
accelerometer rows. In the limit, the nominal accelerometer attitude rows are
field-axis compatible (zero-BA collinearity).

This does **not** contradict the existing jerk-limited collinear-cadence
lemma. That lemma constrains physical specific force. PER-2 constrains the
estimator nominal force appearing in J_att. The proof ledger already contains
an admitted history refuting source-uniform pointwise AW tracking, so no
existing theorem transfers the physical jerk/cadence floor to the nominal
force rows.

Nor does zero homogeneous loss imply zero actual accelerometer innovation:
the compatibility vector is a linearized error direction; the nominal
execution may have nonzero innovations and gains while that homogeneous
direction remains data-null. Therefore the innovation identity cannot turn
PER-2 into physical force collinearity.

Consequently the current assumptions do not rule out persistence of an exact
compatibility line across arbitrarily many words. Conversely, PER-1 does not
construct such an execution; it gives a necessary limiting condition on any
persistent execution.

Thus the persistence question remains undecidable from the current proved
physical-to-nominal relations. What is proved is the dichotomy:
(a) if exact compatibility persists indefinitely, q_W decays geometrically
and the physical tilt scalar has unit persistence, so O2 cannot obtain strict
kernel contraction from BA decay; (b) if the current assumptions exclude the
limiting nominal zero-BA collinearity, the exclusion must come from a new
closed-loop nominal-force reachability theorem, not from MARINE MOTION's
physical jerk lemma alone.

This identifies the same missing bridge as the earlier AW analysis in its
minimal asymptotic form: a relation between physical excitation and the
nominal accelerometer Jacobian on a persistent zero-loss direction.
## Rank-continuous complete-word information representation

Do not form the nuisance Schur complement with a Moore--Penrose inverse.
Retain one fixed-dimensional whitened factor space for each closed event
stratum. Write the complete frozen observation/action record as

`y = O_s x_s + B z`,                                       (RC-1)

where z stacks the nuisance root and every whitened fresh source/noise
factor, and B contains their literal chronological coefficients. Give z its
Euclidean action `|z|^2` by including identity rows for source factors/root
priors that carry finite action. The reduced slow-root quadratic form is

`Q_W(x_s)=min_z ( |O_s x_s+B z|^2 + |R z|^2 )`,             (RC-2)

where R contains exactly the nuisance/source action rows. Equivalently stack

`C_W=[O_s, B; 0, R]`

and minimize the squared norm of `C_W[x_s;z]` over z.

The crucial point is that the nuisance block
`D_W=[B;R]` has a **uniform positive singular floor on its penalized
directions** from the proved nuisance/process/noise floors; exact deterministic
zero-action nuisance directions are retained explicitly in the slow/root
class rather than hidden behind a pseudoinverse. Thus the normal matrix
`D_W'D_W` is positive definite on the chosen nuisance coordinate block after
removing structural zero columns once per event stratum. Its inverse is
continuous. The reduced matrix can therefore be written

`J_W = A_W' [I-D_W(D_W'D_W)^-1 D_W'] A_W`,                  (RC-3)

with `A_W=[O_s;0]`, using an ordinary inverse on the fixed penalized nuisance
space. This representation is continuous through rank changes of the old
observation-only nuisance matrix because those rank changes no longer alter
the rank of D_W.

If a structural deterministic nuisance direction has zero action, move it
into x_s before applying RC-3. The number of such structural coordinates is
fixed by the shipping state architecture, not by the numerical word. Hence
no data-dependent rank stratification or pseudoinverse is needed.

### Direct second-eigenvalue argument

For symmetric PSD J_W,

`lambda_2(J_W)=min_(dim L=2) max_(x in L,|x|=1) Q_W(x)`

and equivalently by Courant--Fischer it is the least information after
allowing one exceptional direction. RC-2 makes Q_W jointly continuous in W
and x on the compact admissible class.

Suppose, for contradiction, `inf_W lambda_2(J_W)=0`. Compactness gives a
subsequence W_n->W_* and, by the min--max characterization, orthonormal
vectors u_n,v_n spanning two-dimensional subspaces with
`Q_Wn(u_n)->0` and `Q_Wn(v_n)->0` (choose eigenvectors for the two smallest
eigenvalues). Passing to subsequences gives orthonormal u_*,v_*.
Continuity of RC-2 yields

`Q_W*(u_*)=Q_W*(v_*)=0`.

Thus `Null(J_W*)` has dimension at least two. But the complete-word zero-action
theorem, including four-S survival and closed MAGNETIC SERVICE on every
boundary event stratum, proves

`dim Null(J_W*)<=1`.

Contradiction. Therefore

`lambda2_bar:=inf_(W admissible) lambda_2(J_W)>0`.           (RC-4)

This proof is continuous through nuisance-observation rank changes and through
appearance/disappearance of the one-dimensional physical compatibility
kernel. It uses only compactness, continuity of the fixed-factor action, and
the already proved nullity<=1 classification.

RC-4 is an existence result; it does not give a numerical lambda2_bar. If any
of the claimed uniform nuisance/action floors needed to make D_W fixed-rank
fails, RC-3 must be replaced by RC-2 directly. The contradiction proof still
works provided minimizers z_n are uniformly bounded modulo structural
zero-action directions; the existing nuisance covariance/action bounds are
the required coercivity input.
## Scalar storage on the persistent compatibility line

Normalize a nontrivial compatibility generator by unit physical attitude:
`nu_W=(theta_hat_W,-q_W)`, `|theta_hat_W|=1`, suppressing zero blocks. A
persistent homogeneous line satisfies

`theta_hat_W -> s_W theta_hat_(W+1)`, `|s_W|=1`,
`q_W -> Phi_b,W q_W = s_W q_(W+1)`.

Consider any positive quadratic scalar storage restricted to this line,
`V_k=lambda^2 w_W`, where

`w_W = a_W + q_W' B_W q_W + 2 c_W' q_W`

is the value of a uniformly positive-definite attitude/BA metric on nu_W.
Uniform equivalence to a fixed physical scalar ceiling requires
`0<w_min<=w_W<=w_max<infinity` on the compact word class.

Along exact persistent transport the scalar amplitude lambda is unchanged
because |s_W|=1. Therefore

`V_(W+1)/V_W = w_(W+1)/w_W`.                                (SI-1)

A source-uniform strict contraction `V_(W+1)<=rho V_W`, rho<1, over an
indefinitely persistent sequence would imply

`w_(W+n)<=rho^n w_W ->0`,

contradicting the required uniform lower equivalence `w_(W+n)>=w_min>0`.

Hence **no uniformly equivalent positive scalar storage can strictly contract
an indefinitely persistent exact compatibility mode with unit physical tilt
amplitude**. Word-dependent BA weighting cannot solve O2 by itself; BA decay
only changes q_W and therefore the representation of the line.

This is a general obstruction, independent of choosing a quadratic metric:
any scalar storage uniformly equivalent to lambda^2 has the same telescoping
contradiction under unit-amplitude persistence.

Therefore O2 can close only if one of the following is proved:
(i) exact compatibility cannot persist indefinitely, so strict-loss words
occur with a source-uniform recurrence; or
(ii) the final practical-stability theorem treats the persistent scalar as a
neutral bounded gauge/class rather than requiring contraction of it. The
latter changes the theorem/state notion and is outside the current declared
regional practical-stability target unless explicitly adopted.

The controlling path under the current theorem is thus persistence exclusion.
A sufficient statement is: there exist N_K<infinity and eta_K>0 such that
within every N_K consecutive MOVING words, at least one word has no exact
compatibility kernel or maps the incoming compatibility line with physical
tilt amplitude <=1-eta_K. This recurring-loss condition would replace
per-word scalar contraction and yield blockwise O2 contraction.

The existing BA decay relation shows that if no such loss occurs, q_W->0.
Thus it suffices to exclude an arbitrarily long sequence whose nominal
accelerometer Jacobians approach the zero-BA field-axis compatibility at all
applied rows. This is exactly the minimal physical-to-nominal excitation
bridge still missing from the proof.
## Infinite persistence and nominal zero-BA collinearity

Assume, for contradiction, an infinite sequence of exact-compatible regular
MOVING words with unit physical tilt persistence. Then q_W->0 geometrically.
At every applied accelerometer epoch the compatibility equation therefore
forces the transverse nominal specific-force Jacobian action toward zero.

The literal mean chronology does not turn this into a contradiction. During
prediction the homogeneous LIN/AW mean follows the OU chain and the S=0
pseudo-update is applied when due **before** the accelerometer correction.
The S update feeds back the estimator's own integrated state; it imposes no
physical measurement constraint. The subsequent accelerometer innovation is
free to correct the AW mean according to the actual measurement residual.

Consequently an asymptotically field-collinear nominal AW history can, at the
level of the mean recursion, be maintained by a sequence of accelerometer
innovations while the S pseudo-feedback removes the resulting integrated
estimator drift. Bounded estimator v,p,S therefore does not imply bounded
physical velocity/displacement or a contradiction with MARINE MOTION.

The exact innovation identity
`f_meas-bhat_a-r_acc=fhat`
also does not help: for prescribed nominal fhat the residual r_acc is the
difference to the physical measurement. Current assumptions bound physical
motion/bias but do not impose a source-uniform bound forcing r_acc->0 or
forcing the estimator nominal AW to track physical acceleration pointwise.
That premise was already refuted in the research ledger.

Thus the infinite-horizon OU+S mean recursion supplies no autonomous
obstruction to nominal zero-BA field-axis compatibility. The estimator can
in principle absorb the physical/nominal mismatch in innovations and its
closed-loop mean corrections; the present proof has no theorem preventing
this indefinitely.

Therefore exact compatibility **cannot currently be proved to break** under
MARINE MOTION + IMU BIAS + MAGNETIC SERVICE alone. The desired blockwise
constants N_K,eta_K do not follow from the established assumptions.

More strongly, any proof of recurring breakage must add a premise/theorem
linking physical excitation to the nominal accelerometer Jacobian or to the
long-run innovation supply (for example a cumulative innovation/action bound
that prevents persistent nominal/physical separation). No such bound is
currently part of the admissibility assumptions or proved from shipping.

This is a theorem-level obstruction, not a missing compactness argument:
compactness can upgrade pointwise eventual breakage to uniform blockwise
breakage only after eventual breakage itself is established. The current
closed-loop equations do not establish it.
## Variational kernel-augmented complete-word coercivity

Fix a retained radius/scalar pair (r,c), c>0. Use the rank-continuous
fixed-factor complete-word representation

`Q_red,W(x)=min_z ||A_W x+D_W z||^2`,

with every literal nuisance/root/source action retained in the fixed factor
space. Let nu_W denote a normalized generator of the (possibly limiting)
one-dimensional physical tilt/BA compatibility relation; when the literal
data nullspace is trivial, retain the closed compatibility-line relation from
the same magnetic/accelerometer equations rather than redefining nu by a
discontinuous least-eigenvector choice.

**Lemma VC.** On the compact admissible complete-word class W(r,c),

`g_*(c,r):=inf_(W, |x|=1)
 [Q_red,W(x)+(1/c)(nu_W' x)^2] >0`.                         (VC-1)

**Proof by contradiction.** Suppose the infimum is zero. Choose
`W_n,x_n`, |x_n|=1, and nuisance minimizers z_n with augmented cost ->0.

1. By the finite closed event/rank stratification and the Arzela--Ascoli
physical-trace result, the literal coefficient histories are compact. Pass to
a subsequence `W_n->W_*` within one closed admissible boundary stratum.
MAGNETIC SERVICE remains valid because it is a closed inequality formed from
actually applied magnetic events; scheduled S and accelerometer rows survive
in real arithmetic by their positive R floors.

2. Nuisance/process coercivity bounds z_n modulo structural deterministic
zero-action directions. Those structural directions have already been moved
into the slow/root coordinate x. Hence z_n is bounded. Pass to
`z_n->z_*`. Also the unit sphere is compact, so `x_n->x_*` with |x_*|=1.

3. Continuity of the literal fixed-factor matrices gives
`A_Wn x_n+D_Wn z_n -> A_W* x_*+D_W* z_*`. Zero limiting cost therefore
gives literal zero complete-word observation/source action at W_*.

4. The complete-word zero-action/nullspace theorem applies on every closed
admissible boundary stratum: four scheduled S rows kill LIN/AW; closed
MAGNETIC SERVICE plus chronological transport kills the gyro quotient; all
accelerometer rows plus the single BA root leave at most the physical
tilt/BA compatibility line. Therefore

`x_* in span(nu_*)`.                                        (VC-2)

5. The kernel-penalty part of the cost also tends to zero. The compatibility
relation is defined by continuous literal magnetic/accelerometer equations;
its graph is closed. Normalize its nontrivial generator by unit physical
attitude amplitude with a fixed sign convention against the transported
field. Along a subsequence the generators converge to a generator nu_* of
the limiting relation. Hence

`nu_*' x_*=0`.                                               (VC-3)

6. From VC-2, `x_*=a nu_*`. With the declared nonzero normalization of nu_*,
VC-3 gives a=0, so x_*=0, contradicting |x_*|=1.

Therefore VC-1 holds.

**Boundary qualification.** If the limiting compatibility relation is
trivial, step 4 already gives x_*=0 and the contradiction is immediate; no
nu_* is needed. If a boundary event pattern violated the nullity theorem,
that would be the exact failure point. The earlier S/accelerometer SPD and
closed MAGNETIC SERVICE arguments exclude the identified event-loss cases.

Thus the variational augmented information has a source-uniform positive
existence floor `g_*(c,r)>0`. This is nonconstructive: it does not yet supply
a numerical value. It supersedes the older unproved DI-6 assertion only at
the existence level and does not revive LE-5--LE-8.
## Diameter and scalar return from the single variational floor

Let `J_aug(W,c)` be the symmetric matrix representing
`Q_red,W(x)+(1/c)(nu_W'x)^2`. Lemma VC gives

`J_aug(W,c)>=g_*(c,r) I` uniformly.                          (VK-1)

Use the same fixed-factor joint model for the terminal state. After optimal
nuisance/source elimination, write the conditional root-to-terminal map as
`Phi_tilde_W` and the known-root terminal covariance as `Pi_W`. The exact
fictitious-kernel covariance is

`P_nu,W=Pi_W+Phi_tilde_W J_aug(W,c)^-1 Phi_tilde_W'`.       (VK-2)

Define the compact terminal-action constant

`L_T(r):=sup_W ||Pi_W^-1/2 Phi_tilde_W||_2^2 < infinity`.   (VK-3)

Finiteness follows from compactness/continuity of the literal fixed-factor
word and the known-root process floor for Pi. No AG leverage estimate is
used.

Then VK-1--VK-3 give

`P_nu,W <= [1+L_T(r)/g_*(c,r)] Pi_W`.

Therefore a usable existence-level diameter is

`K(c,r):=1+L_T(r)/g_*(c,r) < infinity`.                     (VK-4)

Corollary K consequently gives the recurring linear bound
`rho_0(c,r)<=1-1/K(c,r)` whenever the scalar kernel ceiling c is invariant.

### Exact linked scalar return

Do not multiply K by a separately maximized terminal kernel variance. For the
next-root physical compatibility functional nu_+, define directly

`D_W(c):=nu_+' P_nu,W nu_+`

`=nu_+'Pi_W nu_+
 +nu_+'Phi_tilde_W J_aug(W,c)^-1 Phi_tilde_W'nu_+`.         (VK-5)

and

`D(c,r):=sup_W D_W(c)`.                                     (VK-6)

Compactness plus VK-1 makes D finite and attained. The scalar invariant is
exactly

`D(c,r)<=c`.                                                (VK-7)

A coarse usable upper bound follows from the same variational floor:

`D(c,r) <= d0(r)+L_nu(r)/g_*(c,r)`,                         (VK-8)

where
`d0(r)=sup_W nu_+'Pi_W nu_+` and
`L_nu(r)=sup_W ||Phi_tilde_W'nu_+||^2`
with the physical kernel-functional normalization fixed consistently.
Both constants are finite by compactness.

Existence of positive g_* alone does **not** imply VK-7 for some c. As c
increases, the kernel penalty 1/c weakens, so g_*(c,r) can decrease and the
second term can grow proportionally to c on an exact persistent kernel.
The earlier physical-persistence analysis shows that this is a real possible
obstruction, not a defect of VK-8.

Thus O1 is now closed at the existence level by VK-4. O2 is exactly the
one-dimensional fixed-point problem VK-7. The next mathematical task is to
study `D(c,r)/c` using the physical compatibility-line transfer, not to
improve the transverse information floor. If exact unit-persistent kernels
are admissible, `limsup_(c->infinity)D(c,r)/c>=1`; if every admissible chain
breaks compatibility with recurring loss, a strict fixed point may exist.
## Small/large-c analysis of the exact scalar return

For one word W write
`d_W=nu_+'Pi_W nu_+`, `b_W=Phi_tilde_W'nu_+`, and
`J_c=J_W+(1/c)nu_W nu_W'`. Then

`D_W(c)=d_W+b_W'J_c^-1 b_W`.                               (DC-1)

Use coordinates `e_1=nu_W/|nu_W|` and split
`J=[[a,h'];[h,B]]`, `b=(beta,g)`. The rank-one prior adds
`|nu_W|^2/c` to a. Since the augmented matrix is positive for every c>0,
block inversion gives

`D_W(c)=d_W+g'B^-1 g
 + [beta-h'B^-1 g]^2 /
   [a-h'B^-1 h+|nu_W|^2/c]`,                               (DC-2)

with the usual shorted interpretation if B is represented variationally.
The complete-word nullity<=1 theorem and transverse floor make B uniformly
positive on the physical complement.

Define
`j_W:=a-h'B^-1h>=0`,
`ell_W:=beta-h'B^-1g`,
`dperp_W:=d_W+g'B^-1g`. Then exactly

`D_W(c)=dperp_W + ell_W^2/(j_W+|nu_W|^2/c)`.               (DC-3)

This exposes the entire c-dependence.

**Small c.** As c->0+,

`D_W(c)->dperp_W>0` in general, so
`D_W(c)/c -> +infinity` whenever dperp_W>0.                 (DC-4)

Thus very strong fictitious kernel precision cannot satisfy the scalar
invariant because the known-root/process return remains finite and positive.

**Large c, observed word.** If `j_W>0`,
`D_W(c)->dperp_W+ell_W^2/j_W<infinity`, hence
`D_W(c)/c->0`.                                              (DC-5)

Such a word eventually satisfies its scalar inequality.

**Large c, exact-kernel word.** If `j_W=0`, then

`D_W(c)=dperp_W+c ell_W^2/|nu_W|^2`,
`D_W(c)/c=ell_W^2/|nu_W|^2+dperp_W/c`.                     (DC-6)

The limiting slope is

`alpha_W:=ell_W^2/|nu_W|^2`.                               (DC-7)

With the physically normalized compatibility functional, alpha_W is exactly
the squared physical scalar transfer discussed in PT-4/PT-5. Therefore a
persistent exact compatibility line has `alpha_W=1`; a strictly lossy line
has alpha_W<1.

For the source-uniform return
`D(c,r)=sup_W D_W(c)`, compactness makes the supremum attained. Equations
DC-3--DC-7 imply

`limsup_(c->infinity) D(c,r)/c
 = sup_(W with j_W=0) alpha_W`,                             (DC-8)

provided the compact word class and coefficients are fixed at radius r;
words with j_W>0 contribute zero slope. This follows uniformly from the
transverse lambda_2 floor and compactness.

Hence there are two cases.

1. If `alpha_bar(r):=sup_(exact-kernel words) alpha_W <1`, then choose
`epsilon=(1-alpha_bar)/2`. Compactness gives C<infinity such that for all
`c>=C`, `D(c,r)/c<=alpha_bar+epsilon<1`. O2 closes for a finite c.

2. If an admissible exact-kernel word has physical unit persistence
`alpha_W=1`, then DC-6 gives
`D_W(c)/c=1+dperp_W/c>1` for every finite c whenever dperp_W>0.
Therefore `D(c,r)>c` for every finite c: the scalar invariant cannot close
by this one-word ceiling. Its infimum ratio tends to one from above as
`c->infinity`.

Thus the exact O2 question is now completely reduced to whether the compact
admissible word class contains an exact-kernel word with alpha_W=1. No
intermediate finite-c crossing can rescue the invariant in that case.
## Does an admissible one-word exact kernel attain alpha=1?

For an exact-kernel word W, the homogeneous data-null trajectory has zero
gyro-bias quotient, zero LIN/AW root after the four-S argument, and attitude
transport

`theta_N=F_W theta_0`.

In the world-frame historical coordinates F_W is norm-preserving on this
surviving attitude mode. With unit physical-attitude normalization at the
root, `|F_W theta_hat_W|=1`.

However alpha_W in the scalar return is a transfer into the **next-root
compatibility functional** nu_+, not merely the norm of terminal attitude.
Thus alpha_W=1 requires the terminal image of the current null line to lie
exactly on a nontrivial next-root compatibility line with the same unit
attitude amplitude. Equivalently, for some q_W,q_+,

`F_W theta_hat_W = s theta_hat_+`, `|s|=1`,
`Phi_b,W q_W = s q_+`,                                      (A1-1)

and the next word's magnetic/accelerometer compatibility equations must hold
for `(theta_hat_+,-q_+)`.

The current word's exact-kernel equations determine q_W from its own
accelerometer chronology. Equation A1-1 determines the candidate q_+ after
BA decay. But admissibility of W alone imposes no equations from the **next
word's** accelerometer chronology. Therefore a single admissible complete
word does not by itself determine whether its terminal image is an exact
kernel of the following word.

If D_W is defined using a supremum over arbitrary admissible next-root
functionals nu_+, then choosing nu_+ to be the normalized terminal image is
not legitimate unless there exists an admissible successor word whose
compatibility line equals that image. Conversely, ruling out alpha=1 requires
a two-word reachability/compatibility theorem.

Hence the question `does the one-word admissible class contain alpha_W=1?`
is ill-posed unless W includes its admissible successor compatibility data.
The correct compact object is an admissible **two-word pair** `(W_0,W_1)`
with shared terminal/root history. Define alpha on pairs by the transfer from
the exact kernel of W_0 into the exact kernel functional of W_1.

On this pair class, alpha=1 is equivalent exactly to persistence of the
compatibility line across the boundary. The previous analysis has neither
constructed nor excluded such a reachable pair under the current assumptions.
Therefore the present proof cannot answer yes or no from one-word data.

This corrects the scalar-return formulation: D(c,r) must take its supremum
over same-history admissible word pairs (or an execution chain), not over a
word and an independently selected next soft/kernel direction. Compactness
extends to the two-word class, but the alpha=1 issue is precisely the
persistent-compatibility reachability question already identified.
## Same-history two-word alpha compactness dichotomy

Let PairClass(r,c) be the compact class of admissible consecutive MOVING word pairs
`(W0,W1)` sharing the literal terminal/root physical, estimator, covariance,
tuner and scheduler history. Restrict to pairs for which W0 has a nontrivial
exact compatibility line; normalize its generator by unit physical attitude.
Define alpha(W0,W1) as the squared physical-attitude amplitude of the W0
kernel image captured by the W1 compatibility line, and set alpha=0 if W1
has no compatible line containing that image.

Literal world-frame attitude transport is norm preserving on the zero-loss
kernel trajectory, while projection onto a unit next compatibility attitude
line cannot increase norm. Therefore

`0<=alpha(W0,W1)<=1`.                                       (TW-1)

Using the closed compatibility relation rather than a discontinuous chosen
eigenvector makes alpha upper-semicontinuous; on the nontrivial-line stratum
it is continuous. Compactness therefore attains the supremum.

Equality alpha=1 holds iff all of the following are true:

1. W0 has an exact compatibility generator
   `nu0=(theta0,-q0)`, `|theta0|=1`;
2. its terminal attitude image is exactly the unit attitude generator of a
   nontrivial W1 compatibility line:
   `F0 theta0=s theta1`, `|s|=1`;
3. shared BA state transport matches that line:
   `Phi_b,0 q0=s q1`;
4. `(theta1,-q1)` satisfies **all** magnetic and accelerometer zero-loss
   equations of W1.

Conditions 1--4 are precisely an exact compatibility line persisting across
the shared word boundary. There is no additional inequality in the current
assumptions that contradicts them. BA decay only enforces
`|q1|<|q0|` (unless q0=0); it does not reduce the unit attitude amplitude.
MAGNETIC SERVICE is compatible with a one-dimensional field-axis null line
and controls the transverse/gyro sector. Four-S controls LIN/AW. MARINE
MOTION constrains physical attitude span, not the nominal-force compatibility
equations defining q0,q1.

Therefore the existing theorem set cannot prove `sup alpha<1`: excluding the
equality set would require the still-missing physical-to-nominal/reachability
bridge. Conversely, compactness alone does not prove `sup alpha=1`; it only
attains the supremum whatever its value is. An equality pair must still be
shown to exist.

Thus the sharp dichotomy is mathematically well posed but unresolved by the
current assumptions/lemmas:

`alpha_bar:=max_(PairClass) alpha` exists and lies in [0,1];
`alpha_bar<1` iff no admissible exact persistent pair exists;
`alpha_bar=1` iff an admissible exact persistent pair exists.

This is an equivalence, not yet a proof of either branch. It prevents a
compactness overclaim: attainment does not determine whether the attained
maximum equals one.
## Same-history two-word persistence boundary-value problem

Assume conditions 1--4 of TW hold. Pull every accelerometer compatibility
equation in W0 and W1 to the shared boundary world frame. For applied epoch k
in word j, write

`C_jk theta_j = B_jk q_j`,                                  (SH-1)

where C_jk is the literal nominal-force attitude row after chronological
transport and B_jk is the BA transport/row. Persistence gives
`theta_1=F_0 theta_0` with unit physical amplitude and
`q_1=Phi_b,0 q_0` up to the signed frame transport.

Thus the entire two-word zero-loss system is a homogeneous finite-dimensional
boundary-value problem for the single pair (theta_0,q_0). BA decay changes
the right-hand amplitude between words, while the left-hand nominal-force
rows are generated by the actual closed-loop mean history.

Subtracting or comparing SH-1 across the boundary does **not** produce a
contradiction from the current assumptions. The rows C_1k are not prescribed
continuations of C_0k: accelerometer innovations, S feedback, attitude
injection/reset and OU prediction change the nominal AW/attitude history.
The same-history condition couples these rows through shipping dynamics, but
the physical admissibility assumptions do not bound the cumulative innovation
strongly enough to prevent the nominal rows from adapting to the decayed q1.

More explicitly, at each accepted accelerometer correction the physical
measurement supplies a 3-vector innovation. The compatibility constraint
imposes two transverse scalar conditions on the nominal AW/attitude geometry.
Where the transverse AW gain has rank two these conditions are locally
solvable; where it does not, prior reachability work has neither proved nor
excluded the required state. The remaining longitudinal innovation and
endogenous S feedback evolve the internal chain. Therefore the two-word
boundary equations are not overdetermined by a simple dimension count.

MAGNETIC SERVICE likewise does not exclude SH-1: its two-dimensional
information floor is transverse to the allowed field-axis line. MARINE
MOTION controls the physical trajectory but, absent a cumulative
physical-to-nominal innovation bound, does not constrain C_jk enough to make
the BA decay mismatch impossible.

**Theorem status.** Under the current assumptions and proved shipping
invariants, conditions 1--4 are neither shown reachable nor excluded. The
same-history formulation removes independent-row freedom but does not add a
proved inequality that forces inconsistency. Therefore no rigorous theorem
`alpha_bar<1` or `alpha_bar=1` follows.

A genuine exclusion theorem would require one additional proved consequence
of existing shipping/admissibility, for example a uniform cumulative bound
on the transverse accelerometer innovation/action over two words that is
strictly smaller than the amount required to retune compatibility from q0 to
`Phi_b q0`; or a closed-loop reachability invariant forbidding that retuning.
No such consequence is currently established.

Conversely, a reachability theorem requires solving the literal two-word mean
boundary-value system and verifying MARINE MOTION and strict MAGNETIC SERVICE;
the local AW-rank argument alone is insufficient.

This is the terminal reachability/exclusion gap for O2. Further covariance
or compactness arguments cannot decide it because both branches are
compatible with all presently proved inequalities.
## Controlling O1 path: complete-word nullspace and closed-range continuity

All O1 work is consolidated to the path
`complete-word nullspace -> uniform qualitative augmented floor -> quantitative enclosure`.
The later detectability/AW-gain/sign-crossing/soft-kernel investigations are
diagnostic or alternative formulations and are not controlling.

### Theorem A: complete-word slow nullspace

For every admissible recurring complete A21 word W, including every closed
event-boundary word allowed by actually-applied MAGNETIC SERVICE,

`Ker J_s(W) subset span(nu_W)`.                              (A)

Here J_s is the nuisance-eliminated slow-root information defined
variationally, and nu_W is the physical tilt/BA compatibility generator when
that relation is nontrivial; if the compatibility relation is trivial,
`Ker J_s(W)={0}`.

The proof is the literal zero-action classification already established:
zero nuisance/process action rigidifies one deterministic trajectory; four
scheduled S rows kill LIN/AW; magnetic service plus chronological transport
eliminates gyro-bias quotient directions; all accelerometer rows plus the
single BA root leave at most the physical tilt/BA compatibility line.
Scheduled S and accelerometer rows survive boundary limits in real arithmetic
because their innovation covariances have positive R floors; MAGNETIC SERVICE
is closed and uses actually applied rows. Thus the same classification holds
on admissible boundary words.

### Closed-range convergence lemma

Let
`F_W=Sigma_W^-1/2 O_f,W`
be the whitened nuisance observation operator in a fixed event stratum.
The desired graph property

`W_n->W, z_n in Range(F_Wn), z_n->z => z in Range(F_W)`     (CR-1)

is **not true for arbitrary continuous matrices** when rank drops: for
`F_epsilon=diag(1,epsilon)`, the ranges are R^2 for epsilon!=0 but limit to
`span(e1)` at epsilon=0, while z_n=e2 violates CR-1.

Therefore observation-range closedness cannot be assumed across nuisance-rank
loss. This is exactly why a pseudoinverse/projector representation is unsafe.

For the literal proof, avoid CR-1 by retaining nuisance/process action in the
fixed-factor variational form

`Q_red,W(x)=min_y ||A_W x+D_W y||^2`,                        (CR-2)

where D_W stacks both nuisance observation columns and their process/root
action rows. Existing nuisance/process coercivity gives a uniform lower
singular floor for D_W after structural zero-action directions are moved into
the slow root. Hence Range(D_W) has fixed rank and a closed graph, minimizers
are uniformly bounded, and Q_red is continuous (in particular lower
semicontinuous) on each closed stratum and across the admissible finite union.

If a claimed nuisance coordinate lacks such an action floor, it must be moved
into the slow/root state; otherwise the compactness proof fails precisely at
that structural rank-loss direction. No observation-only closed-range lemma
is used.

### Uniform qualitative O1

Fix `mu=1/c>0`. Suppose no uniform floor exists. Then there are admissible
`W_n` and `|x_n|=1` with

`x_n'[J_s(W_n)+mu nu_n nu_n']x_n ->0`.                      (UO-1)

Compactness gives, after subsequences, `W_n->W_*`, `x_n->x_*`,
`|x_*|=1`. The fixed-factor lower semicontinuity from CR-2 gives

`x_*'J_s(W_*)x_*=0`.                                        (UO-2)

By Theorem A, `x_*=alpha nu_*` (or x_*=0 immediately if the limiting
compatibility relation is trivial). The kernel penalty in UO-1 tends to zero.
The physical compatibility relation has a closed graph; with its fixed
normalization, pass to a subsequence `nu_n->nu_*`, obtaining

`nu_*'x_*=0`.                                                (UO-3)

Thus `alpha=0`, so x_*=0, contradicting |x_*|=1. Therefore

`J_s(W)+mu nu_W nu_W' >= g_mu(r,c) I`, `g_mu(r,c)>0`,       (UO-4)

uniformly over all admissible complete A21 words.

This is qualitative only. The next obligation is a constructive enclosure
`g_mu>=g_under(c,r)>0`; only after that should K(c,r), the same-history scalar
return D(c,r), rho_0, and the nonlinear retained radius be developed.
## Constructive enclosure for the variational augmented floor

Quantify Theorem A's implication chain in the same fixed-factor action. Let
`E` denote total variational action including the kernel row. Introduce
sector amplitudes u=(u_N,u_L,u_G,u_K), where u_N is penalized nuisance/source
amplitude, u_L the homogeneous LIN/AW root amplitude, u_G the magnetic/gyro
quotient amplitude, and u_K the physical tilt/BA compatibility amplitude.
Use the literal same-history maps; these amplitudes are bookkeeping
coordinates, not independent worst-case boxes.

1. **Nuisance/process coercivity.** Existing process/noise floors give
`E >= a_N^2 u_N^2`, with constructive `a_N>0` from the minimum whitened
fresh-factor/root-action floor on the retained class.

2. **Four-S LIN/AW implication.** Four scheduled S rows give
`|R_S4 x| >= a_L u_L - b_L u_N`.                            (CG-1)
The qualitative Chebyshev argument makes a_L positive; constructively one
may use the generalized-Vandermonde determinant with the correct tau_min
dependence plus explicit row-norm ceilings, or interval-enclose the analytic
minimum. The previous false >6.82e4 bound is not used.

3. **Magnetic/gyro implication.** MAGNETIC SERVICE + Lemma T give
`|R_MG x| >= a_G u_G - b_G,L u_L - b_G,N u_N`,              (CG-2)
where a_G is built from mu_M and the chronological transport modulus.

4. **Accelerometer tilt/BA implication.** After the preceding sectors are
charged, all literal accelerometer rows give
`|R_A x| >= a_K dist(x_K,span nu_W)
              -b_K,G u_G-b_K,L u_L-b_K,N u_N`.             (CG-3)
The constructive coefficient a_K is the only genuinely new quantitative
modulus: it is the minimum nonzero singular value of the complete
same-history magnetic-compatible accelerometer/BA map on the compact
retained class, **after allowing the physical line nu_W**. It is not the
retired absolute principal-angle floor; the line is retained explicitly.

5. **Kernel row.** With mu=1/c,
`sqrt(E) >= sqrt(mu) |nu_W' x|`.                            (CG-4)
Together with CG-3 this controls the full tilt/BA sector.

Collect the inequalities before applying triangle/Cauchy into the lower
triangular implication matrix

`M_imp = [[a_N,0,0,0],
          [-b_L,a_L,0,0],
          [-b_GN,-b_GL,a_G,0],
          [-b_KN,-b_KL,-b_KG,a_K]]`,

and append the kernel functional row `sqrt(mu) k_W'` on the final physical
sector. Let `H_imp(c,r,W)` be the resulting small Gram matrix. Then

`E >= u' H_imp u`,
`g_under(c,r):= inf_W lambda_min(H_imp(c,r,W))`.             (CG-5)

If the coefficients have source-uniform bounds with positive diagonal
`a_N,a_L,a_G,a_K`, compactness gives an explicit computable
`g_under>0`. This is a constructive enclosure of g_mu and uses no
observation pseudoinverse or unprojected gyro leverage.

**Current quantitative status.** a_N and the magnetic-service part of a_G
come from existing explicit floors. a_L is analytically extractable from the
correct four-S determinant/norm bounds (with tau_min, not tau_max). The
remaining coefficient a_K has not yet been numerically/analytically
lower-bounded in the current proof. Qualitative Theorem A proves it cannot
vanish in the full augmented implication chain, but a separate positive
minimum for CG-3 alone may fail near rotation of the physical compatibility
line. Therefore CG-5 must, if necessary, keep the combined accelerometer +
kernel rows as one block rather than require a_K>0 separately.

The safe constructive target is thus the minimum singular value of that
**combined literal block** on the compact class. It is an analytic
finite-dimensional continuous minimization with the physical line retained;
extracting a certified lower enclosure for it is the remaining numerical-
constant task before K(c,r) can be evaluated.
## Combined accelerometer/kernel block floor

After charging nuisance, LIN/AW and magnetic/gyro sectors, let A_W denote the
literal reduced accelerometer map on the remaining tilt/BA slow coordinates.
Theorem A gives

`Ker A_W subset span(nu_W)`

on this reduced sector; when the physical compatibility line is nontrivial,
equality is the worst case. Define

`B_W(c)=[A_W; c^-1/2 nu_W']`.

Then
`B_W'B_W=A_W'A_W+c^-1 nu_W nu_W'`.                          (AK-1)

For every word this matrix is positive definite: if B_W x=0, then
`A_Wx=0`, so x=a nu_W, while `nu_W'x=0` forces a=0.

To obtain a quantitative formula, normalize `e=nu_W/|nu_W|` and decompose
`x=t e+y`, y perpendicular e. Write the reduced accelerometer Gram as

`A_W'A_W=[[a,h'];[h,C]]`

in `(e,e_perp)` coordinates. Since its nullspace is contained in span(e),
`C>0`. The augmented Gram is

`G_AK=[[a+|nu|^2/c,h'];[h,C]]`.                             (AK-2)

Its Schur complement is

`s_AK = a+|nu|^2/c-h'C^-1h
     = j_A(W)+|nu|^2/c`,                                   (AK-3)

where `j_A>=0` is the accelerometer information shorted to the physical line.

Let
`gamma_perp(W)=lambda_min(C_W)>0` and
`H_A(W)=||C_W^-1/2 h_W||`. A direct completion of squares gives, for
`x=(t,y)`,

`x'G_AK x
 = (y+C^-1 h t)' C (y+C^-1 h t) + s_AK t^2`.               (AK-4)

Using `|y+C^-1ht|^2 + t^2` versus `|y|^2+t^2` yields the explicit bound

`lambda_min(G_AK)
 >= min(gamma_perp,s_AK) / (1+||C^-1 h||)^2`.               (AK-5)

(a sharper 2x2 norm-equivalence factor may be substituted). Therefore

`sigma_min(B_W(c))
 >= sqrt(min(gamma_perp(W), j_A(W)+|nu_W|^2/c))
    /(1+||C_W^-1 h_W||).                                    (AK-6)

Now use compactness. The physical generator normalization gives
`|nu_W|>=nu_min>0`; the reduced transverse accelerometer block C_W is
continuous and has no zero direction by Theorem A after the previously
charged sectors, so
`gamma_bar_perp=inf_W gamma_perp(W)>0`; and continuity gives
`Hbar_A=sup_W||C_W^-1 h_W||<infinity`. Hence

`a_AK_under(c,r)
 := sqrt(min(gamma_bar_perp,nu_min^2/c))/(1+Hbar_A) >0`.     (AK-7)

This is a constructive symbolic enclosure:
`a_AK_under <= inf_W sigma_min B_W(c)`.

Important qualification: the assertion `gamma_bar_perp>0` is valid only if
the reduced A_W used here is the **complete same-history accelerometer block
after** nuisance/LIN/gyro implications have been incorporated, so that any
additional zero direction would contradict Theorem A. It must not be
identified with an isolated two-epoch or pointwise accelerometer matrix; that
would resurrect the refuted principal-angle route.

Thus the combined block has a certified positive enclosure in terms of three
compactness constants `(gamma_bar_perp,nu_min,Hbar_A)`. Numerical/closed-form
values for those constants are still to be extracted from the literal
coefficient bounds before H_imp and K can be evaluated numerically.
## Explicit implication-matrix floor and Riccati diameter

Use the sector implication operator

`M=[[aN,0,0,0],
    [-bLN,aL,0,0],
    [-bGN,-bGL,aG,0],
    [-bKN,-bKL,-bKG,aAK,]]`,                               (IG-1)

where `aAK=a_AK_under(c,r)` from AK-7. All coefficients are nonnegative
source-uniform enclosures on the retained class; the b's bound the linked
same-history leakage terms from already charged sectors.

The total fixed-factor action satisfies
`E >= |M u|^2` after the sector amplitudes are normalized consistently.
Therefore

`g_under(c,r)=sigma_min(M)^2`                               (IG-2)

is a valid constructive augmented-information floor.

To make positivity explicit without a diagonally-dominant assumption, bound
the inverse triangular matrix. Let

`d1=1/aN`,
`d2=(1+bLN*d1)/aL`,
`d3=(1+bGN*d1+bGL*d2)/aG`,
`d4=(1+bKN*d1+bKL*d2+bKG*d3)/aAK`.                         (IG-3)

These recursively bound the row-sum norm of the solution of `M u=y`:
`|u_i|<=d_i ||y||_infinity`. Hence

`||M^-1||_2 <= 2 sqrt(d1^2+d2^2+d3^2+d4^2)`                (IG-4)

(the factor 2 converts the four-dimensional infinity norm of y to its
Euclidean norm; a sharper norm conversion may be used). Consequently

`sigma_min(M) >=
 1/[2 sqrt(d1^2+d2^2+d3^2+d4^2)]`,                         (IG-5)

and the explicit symbolic floor is

`g_under(c,r):=
 1/[4(d1^2+d2^2+d3^2+d4^2)] >0`.                           (IG-6)

Because
`aAK(c,r)=sqrt(min(gamma_bar_perp,nu_min^2/c))/(1+Hbar_A)`,
all c-dependence is explicit in d4. No Gershgorin subtraction or
cross-coupling smallness assumption is needed.

The variational theorem gives `J_aug>=g_under I`. With
`L_T(r)=sup_W ||Pi_W^-1/2 Phi_tilde_W||^2`, the exact terminal covariance
identity yields

`P_nu,W <= K(c,r) Pi_W`,

`K(c,r):=1+L_T(r)/g_under(c,r)`.                            (IG-7)

Hence

`rho_0(c,r)<=1-1/K(c,r)<1`

whenever the same-history scalar kernel ceiling is invariant.

All constants in IG-3 are now named proof obligations rather than hidden
matrix minima. Existing results supply aN and the magnetic-service component
of aG; aL requires the corrected four-S quantitative enclosure; aAK is AK-7;
the b coefficients are literal operator-norm bounds of the linked leakage
maps and must be extracted from the retained coefficient bounds. Until those
numbers/closed forms are supplied, IG-6 is an explicit symbolic theorem, not
a numerical certificate.

The controlling next step is D(c,r)<=c on admissible same-history word pairs.
No additional O1 architecture is needed.
## Constant ledger for the fixed O1 architecture

The implication architecture is fixed. The constants split into three classes.

**Already explicit from shipping/proved floors.**
`a_N` is the minimum penalized nuisance/process whitening floor after
structural zero-action coordinates are moved to the slow root. `a_G` has an
explicit magnetic-service factor from `mu_M` multiplied by the proved
chronological-transport/Lemma-T modulus. `nu_min>=1` under unit physical
attitude normalization of `nu=(theta_hat,-q)`.

**Analytically extractable from bounded literal matrices.**
`a_L` comes from the four-S generalized Vandermonde with the correct
`tau_min` dependence and explicit row-norm ceiling. The six b coefficients
are operator norms of the literal same-history leakage blocks and admit
finite enclosures from the retained dt/state/tuner/reset/process bounds.
`Hbar_A` and `L_T` are likewise finite operator-norm enclosures once the
known-root process floor for Pi is inserted.

**Still only compactness-positive unless further quantified.**
`gamma_bar_perp` is the transverse minimum eigenvalue of the complete reduced
accelerometer block. Theorem A + compactness proves it positive in this
representation, but no closed numerical lower enclosure is currently
derived. Therefore the current `g_under` and K are constructive symbolic
forms but not yet evaluable numerical certificates.

Do not introduce another O1 route to quantify gamma_bar_perp; if a numerical
certificate is required, interval/analytic enclosure of this fixed literal
block is the permitted constant-extraction task.

## Same-history scalar return fixed point

Now consider an admissible shared-history pair `(W0,W1)`. Let `nu0` be W0's
physical compatibility functional and `nu1` W1's. For W0 define

`d_01=nu1' Pi_0 nu1`,
`b_01=Phi_tilde_0' nu1`,
`J0(c)=J_s(W0)+(1/c)nu0 nu0'`.

Then exactly

`D_01(c)=d_01+b_01' J0(c)^-1 b_01`,                         (DP-1)
`D(c,r)=sup_(same-history pairs) D_01(c)`.                  (DP-2)

Split J_s(W0) along the normalized physical line e0=nu0/|nu0|. Shorting the
transverse block gives

`D_01(c)=dperp_01 + ell_01^2/(j_0+|nu0|^2/c)`,              (DP-3)

with `j_0>=0`. Hence D_01 is continuous, nondecreasing and concave in c.

If `j_0>0`, D_01(c) is bounded as c->infinity and eventually lies below c.
If `j_0=0`,

`D_01(c)=dperp_01+c alpha_01`,
`alpha_01=ell_01^2/|nu0|^2`.                               (DP-4)

On the same-history pair class alpha_01 is exactly the physical scalar
transfer into W1's compatibility functional. Therefore:

- if `alpha_bar:=sup_(exact-kernel same-history pairs) alpha_01 <1`,
  compactness and the uniform transverse floor give a finite C such that
  `D(c,r)<c` for all sufficiently large c;
- if an admissible same-history pair has `alpha_01=1` and `dperp_01>0`,
  then `D(c,r)>c` for every finite c.

Thus the existence of a finite scalar fixed point is equivalent to excluding
a unit-persistent exact-kernel pair (or otherwise proving its dperp term
vanishes, which the known-root process covariance generally prevents).

The current reachability work has neither constructed nor excluded such a
pair. Therefore **no finite c is presently certified**, and it would be
premature to extract numerical values for all O1 constants as though O2 were
known to close.

The next controlling obligation remains the same-history equality case
`alpha_01=1`: prove it unreachable under the literal closed-loop mean
dynamics, or construct it. Only the first outcome permits a finite D<=c
certificate in this one-word scalar-ceiling architecture.
## Two-word finite-horizon detectability at the unit-persistence equality case

Focus only on a same-history pair `(W0,W1)` with W0 exact kernel and suppose
`alpha_01=1`. Let x0=lambda nu0 be the W0 zero-action mode with unit physical
tilt normalization. Its deterministic terminal image is x1=lambda nu1 if
equality holds. W0 contributes zero complete-word action by definition.

Define the literal **joint two-word action**

`Q_2(x0)=min_(all shared nuisance/source factors)
 [action on W0 + action on W1]`,                            (FD-1)

with the W0 terminal state constrained to be the W1 root; do not re-minimize
that shared state independently. A source-uniform finite-horizon detectability
inequality sufficient for O2 is

`|lambda|^2 <= C_FD Q_2(x0)`                               (FD-2)

for every unit-persistence candidate pair. If FD-2 holds, Q_2 cannot vanish
on a nonzero equality mode, so alpha_01=1 is impossible for an exact
two-word null trajectory; compactness then gives a strict loss margin.

Now classify zero joint action. Zero W0 action gives x0 in span(nu0). Zero
fresh/process action across the shared boundary propagates the BA component
homogeneously by Phi_b,0 and the attitude component by the norm-preserving
world-frame transport. Zero W1 measurement action then requires exactly

`Phi_b,0 q0 = s q1`, `|s|=1`,                              (FD-3)

together with every W1 accelerometer compatibility equation. Thus **zero
joint action is equivalent to the same-history persistent compatibility
conditions 1--4** already identified.

Therefore finite-horizon detectability FD-2 is not a consequence of BA decay
alone: BA decay is already incorporated in FD-3 and can be absorbed by the
changed next compatibility ratio q1. To prove FD-2 one must show that the
literal W1 nominal-force rows cannot all satisfy the retuned compatibility
equations while sharing the actual W0 terminal/mean history.

The closed-loop mean equations currently supply no such inequality. The
accelerometer innovations that generate the nominal AW/attitude history are
not part of the homogeneous error action Q_2 as penalized disturbances;
they belong to the base trajectory. MARINE MOTION bounds the physical base
trajectory but no proved cumulative innovation bound transfers that physical
excitation to the nominal W1 Jacobian rows. Consequently a persistent
compatibility pair is not excluded by the current joint-action algebra.

**Exact conclusion.** A source-uniform FD inequality of form FD-2 is
equivalent, on the equality class, to excluding conditions 1--4. The current
assumptions/lemmas do not prove it. Recasting the problem as joint action
does not create the missing physical-to-nominal bridge.

A noncircular detectability theorem would need to enlarge the action with a
quantity already bounded by admissibility and sensitive to the base
accelerometer innovation/nominal-force retuning, then prove that maintaining
FD-3 for two words consumes a positive amount of that quantity. No such
base-trajectory supply term is currently present in the stability action.

Thus the equality case remains the sole O2 blocker: neither unreachable nor
constructed. Further homogeneous covariance/information estimates cannot
settle it.
## Base innovation supply versus persistent compatibility

Use the exact accepted-accelerometer identity
`r_acc,k=f_meas,k-bhat_a,k-fhat_k`.                         (BI-1)
Let P_bperp denote projection transverse to the transported physical magnetic
axis. Along a zero-BA limiting compatibility sequence,
`P_bperp fhat_k ->0` at every applied accelerometer epoch. Therefore

`P_bperp r_acc,k
 =P_bperp(f_meas,k-bhat_a,k)+o(1)`.                         (BI-2)

The physical jerk/cadence lemma says dense applied epochs cannot all have
physical specific force parallel to b on a MOVING interval. With the bounded
physical accelerometer bias, this supplies nonzero transverse **innovation**
at some epochs unless the physical transverse force is canceled by bhat_a.
However bhat_a is an estimator mean, not the bounded physical bias; its
shipping projection only bounds its magnitude.

At an accepted correction the nominal AW/attitude means change by
`Delta a_w=K_aw r_acc`, `Delta theta=K_theta r_acc`.        (BI-3)
To force a break of compatibility from BI-2 one needs a lower gain/action
inequality such as

`|P_bperp[K_aw;K_theta] r| >= kappa_I |P_bperp r|`          (BI-4)

on the relevant innovation directions, or a multi-epoch analogue after S
feedback. The prior AW-gain reachability analysis proved no such uniform
lower singular value for reachable A21 covariances; cross-covariances can
algebraically cancel the AW gain block, and reachability of cancellation is
open.

More importantly, even BI-4 would not by itself force incompatibility:
persistent compatibility can use the correction to **retune** the nominal
force to the next decayed BA ratio. The required retuning is O(1-Phi_b) per
word and therefore finite. There is no admissibility bound on cumulative
base innovation energy `sum r_acc' S^-1 r_acc`; NIS is a gate/diagnostic, not
a stated finite long-run budget. Thus persistent retuning need not exhaust a
proved resource.

Consequently the existing base innovation/physical-motion relationship does
not exclude a persistent pair. The physical jerk lemma guarantees excitation
of the measurement, but shipping corrections are designed to absorb such
innovation, and the proof assumptions place no cumulative budget on that
absorption.

### Constructive-pair status

Conversely, BI-1--BI-3 make a persistent pair locally plausible where the
transverse correction gain has rank two: the physical measurement innovation
provides exactly the controls needed to enforce the two transverse
compatibility equations, while BA decay merely changes the target ratio.
But a rigorous pair construction still requires solving the endogenous
longitudinal/S/BA endpoint equations and proving the resulting sampled
physical measurements admit a MARINE MOTION continuation with strict
MAGNETIC SERVICE. Those obligations remain unsolved.

Therefore neither branch closes under the current theorem set. What is now
proved is that no existing physical-motion bound creates a finite innovation
supply budget capable of excluding persistent compatibility. Any exclusion
theorem needs a new consequence of shipping/admissibility beyond BI-1--BI-3;
otherwise the productive direction is an explicit boundary-value construction
of the persistent pair.
## Explicit two-word compatibility boundary-value construction

Fix one regular two-word scheduler/event pattern with strict magnetic-service
margin and no projection/gate boundary. Let the accepted accelerometer
innovations `u_k=r_acc,k in R3` be the base-trajectory controls. The literal
mean recursion over the pair is a smooth finite-dimensional map

`m_(k+1)=F_k(m_k,u_k)`,                                     (BV-1)

where prediction includes endogenous S pseudo-feedback when due; S residuals
are not independent controls.

At each accelerometer epoch impose the two transverse compatibility equations

`C_k(m_k,u_<k) theta_k - B_k(m_k,u_<k) q_k =0`.             (BV-2)

Whenever the transverse derivative with respect to two components of u_k has
rank two, the implicit-function theorem solves
`u_k,perp=psi_k(m_k,u_k,parallel)` exactly. Substitution leaves one scalar
control `ell_k=u_k,parallel` per accelerometer epoch and a reduced compatible
recursion

`m_(k+1)=Fhat_k(m_k,ell_k)`.                                (BV-3)

Across the shared boundary enforce unit-persistent transport
`theta_1=F_0 theta_0`, `|theta_1|=|theta_0|=1`, and
`q_1=Phi_b,0 q_0` in the literal frame. The second word's BV-2 equations are
then imposed by the same elimination.

After transverse elimination there is **no independent BA endpoint closure
condition**: BA mean evolution is part of m and q is the homogeneous error
compatibility ratio, not the estimator BA mean. Likewise S need not return
to its initial value; the two-word pair is not required to be periodic.
Requiring S_N=S_0 or BA_mean,N=BA_mean,0 would add artificial constraints.

Therefore the finite pair construction needs only:
(i) existence of a regular initial A21 state with rank-two transverse control
at each needed epoch (or continuation through isolated rank changes);
(ii) compatibility equations BV-2 through both words; and
(iii) physical realization/service admissibility of the resulting controls.
There is no separate longitudinal/S/BA algebraic endpoint equation forced by
conditions 1--4.

This corrects the earlier boundary-value count: the remaining longitudinal
controls are free parameters, not variables that must close an estimator
periodic orbit. Locally, if a strict-margin compatible base execution with
rank-two transverse gain exists, the IFT propagates an exact-compatible
two-word family without an endpoint overdetermination.

However the construction still lacks that **base execution**. Constructor
diagonal covariance has the required AW gain rank but is not a certified
recurring A21 root; reachable regular A21 rank-two gain was not proved.
Thus the boundary-value equations themselves reveal no longitudinal/S/BA
inconsistency. The only remaining construction obstruction is reachable-base
existence plus physical/service realization.

Consequently no exclusion theorem arises from endpoint closure. If a reachable
strict-margin rank-two A21 base can be exhibited, conditions 1--4 are locally
solvable across two words and alpha_bar=1 follows after physical/service
verification. Without such a base, reachability remains the unresolved step.
## Physical realization of a compatible innovation sequence

Assume a regular strict-margin A21 base execution and the IFT-compatible
innovation sequence u_k from BV-3. The accepted accelerometer identity defines
the required sampled calibrated physical specific force exactly:

`f_phys,k = fhat_k + bhat_a,k + u_k`                         (PR-1)

(with the literal temperature/calibration convention). Choose the physical
sensor bias trace within IMU BIAS and set the physical translational
acceleration samples from PR-1 and the known attitude/gravity transform.

On a finite two-word interval, any sufficiently small perturbation of a
strictly admissible base acceleration sample sequence can be interpolated by
a C1 piecewise-cubic/Hermite acceleration trace matching the samples. Its
jerk norm is bounded by a constant times the maximum sample perturbation
divided by the minimum sample gap. Therefore, if the base execution has
strict margins to the acceleration and jerk limits, the IFT controls can be
restricted to a neighborhood in which the interpolated physical trace
preserves those bounds.

Integrating the acceleration perturbation over the fixed two-word horizon
changes velocity and displacement continuously. Strict base margins to their
bounds therefore persist for sufficiently small controls. The bounded-
potential inequalities on the finite interval are also continuous in the C0
position norm. Outside the two-word interval splice back to the original
globally admissible base continuation with a compact smooth transition;
strict margins and the free longitudinal controls permit zero net velocity/
position perturbation to be imposed if needed for the splice. This gives a
global MARINE MOTION continuation provided the base history has strict
margins.

Keep the physical attitude/gyro/magnetic-field history equal to the base
history when realizing the translational perturbation. The accelerometer
measurement changes, but the physical magnetic samples do not. Estimator
attitude means/resets can nevertheless change because accelerometer
innovations inject attitude. The transported magnetic service rows therefore
vary continuously with the IFT controls rather than remaining exactly fixed.

If the base execution satisfies MAGNETIC SERVICE with strict margin
`lambda_min(G_M)>=mu_M+delta_M`, delta_M>0`, continuity of the finite
transported-row Gram gives an epsilon_M>0 such that all sufficiently small
compatible controls retain
`lambda_min(G_M)>=mu_M`. Thus strict magnetic service is an open property on
a fixed accepted-event branch.

Likewise NIS/finite-input/projection branch conditions with strict margins
are preserved for sufficiently small controls. Therefore **physical
realization and MAGNETIC SERVICE are locally open obligations**, not a new
algebraic obstruction.

The remaining nonlocal issue is the base point: one still needs an actual
reachable recurring A21 execution that simultaneously (i) has strict physical
and service margins, (ii) lies on an exact compatibility line, and (iii) has
rank-two transverse compatibility control. If such a base exists, the IFT
construction plus PR-1 yields a genuine admissible persistent two-word pair
and hence alpha_bar=1. The current proof has not established existence of
that base.
## Fixed-metric exact compatibility-mode return

Fix once and for all a positive-definite physical root metric M on the
attitude/BA slow sector, with declared units/scales. Normalize every nontrivial
compatibility generator by

`nu_W' M nu_W = 1`.                                        (MR-1)

Transport the exact W kernel homogeneously to the next recurring root:

`nu_hat_(W->+) = T_W nu_W`,                                 (MR-2)

where T_W is the literal deterministic attitude/reset + BA OU root transport
on the zero-loss mode. For a same-history successor W+ with normalized
compatibility generator nu_+, define the metric projection coefficient

`a_W = nu_+' M nu_hat_(W->+)`                               (MR-3)

(the denominator is one by MR-1; retain it explicitly for other
normalizations). This coefficient is invariant under coordinate rescaling
when M is transformed accordingly.

The O2 equality question is now the theorem statement

`alpha_M := sup_(admissible exact-kernel same-history pairs) |a_W| < 1 ?`

or, if one-word strict loss fails, whether there exist finite m and delta>0
with every admissible m-word exact-kernel chain satisfying
`|prod a_j|<=1-delta`.

BA decay enters T_W exactly. It does not by itself imply |a_W|<1 because the
next compatibility line may change its BA/attitude ratio. Equality |a_W|=1
means the transported old line is M-collinear with the next compatibility
line with no metric amplitude loss. This is the precise persistent-line
condition to characterize/reachability-test.

Until MR-3 is bounded strictly below one (or a finite product is), no theorem
statement should use Euclidean `|n_+'Phi n|<1`, and no finite scalar O2
invariant is claimed.
## m-word equality-chain theorem for the metric-normalized compatibility mode

Fix the positive physical root metric M and normalize every nontrivial
compatibility generator by `nu_j' M nu_j=1`. For an admissible same-history
chain of complete MOVING words W_0,...,W_(m-1), let T_j be the literal
zero-loss deterministic root transport and define

`a_j = nu_(j+1)' M T_j nu_j`.                               (EC-1)

By M-Cauchy--Schwarz, after normalizing the transported physical attitude
amplitude consistently, `|a_j|<=1`; equality holds exactly when the
transported old compatibility line is M-collinear with the next line with no
metric amplitude loss.

Assume equality at every boundary of a chain. Then there are signs s_j,
`|s_j|=1`, such that

`T_j nu_j = s_j nu_(j+1)`.                                 (EC-2)

Writing `nu_j=(theta_j,-q_j)` in the physical attitude/BA coordinates gives

`F_j theta_j=s_j theta_(j+1)`,
`Phi_b,j R_ba,j q_j=s_j q_(j+1)`.                           (EC-3)

Therefore

`|q_m| <= (prod_(j=0)^(m-1) Phi_b,j) C_R |q_0|`.           (EC-4)

With positive minimum word duration and bounded BA time constant there is a
source-uniform `phi_bar<1` (frame transports are norm preserving in the
physical BA norm), so

`|q_m| <= phi_bar^m |q_0|`.                                (EC-5)

Thus every infinite equality chain satisfies `q_j->0` exponentially. The
closed compatibility equations then imply that every subsequential limiting
word has the zero-BA field-axis compatibility line.

**Compactness consequence.** Let C_infty be the compact inverse-limit class
of admissible infinite recurring MOVING executions (use diagonal
Arzela--Ascoli/event compactness). If no execution in C_infty can satisfy
`|a_j|=1` for all j, then there exist finite m and delta>0 such that every
admissible m-word chain obeys

`|prod_(j=0)^(m-1) a_j| <= 1-delta`.                        (EC-6)

Proof: otherwise for every n choose an n-word chain with product ->1.
Since every factor is <=1, each fixed-prefix factor tends to equality.
Diagonal compactness produces an infinite admissible limiting execution with
`|a_j|=1` at every boundary, contradiction. Continuity/upper-semicontinuity
then upgrades the finite m maximum to a strict `1-delta`.

Conversely, if an infinite equality execution exists, no finite m,delta can
satisfy EC-6.

Hence EC-6 is **equivalent** to excluding an infinite exact persistent
compatibility execution. Equations EC-3--EC-5 characterize every such
execution asymptotically: it must approach the zero-BA field-axis
compatibility manifold while retaining unit M-normalized scalar amplitude.

The current assumptions do not yet exclude that limiting execution. The
literal base innovation analysis shows physical excitation can appear as
accelerometer innovation and can in principle retune the nominal force; no
proved cumulative innovation budget prevents this indefinitely. Therefore
the m-word theorem is proved conditionally:

`no infinite equality chain  <=>  exists finite m,delta with EC-6`,

but the left-hand exclusion is still OPEN. No claim of blockwise O2
contraction is made until that exclusion or an explicit infinite equality
construction is supplied.
## Infinite equality execution: limiting closed-loop analysis

Assume an admissible infinite recurring MOVING execution with `|a_j|=1` at
every word boundary. EC-5 gives `q_j->0` exponentially. Hence at every
applied accelerometer epoch in late words the homogeneous compatibility
equation tends to

`J_att,k theta_k =0`,                                       (IE-1)

where theta_k is the transported physical magnetic-axis tilt. In the usual
accelerometer row this means the estimator nominal specific force tends to
the magnetic-axis-compatible direction.

Use the exact base identity
`r_acc,k=f_meas,k-bhat_a,k-fhat_k`.                         (IE-2)
All terms are bounded on the retained region: physical acceleration/bias are
bounded by admissibility and estimator means are bounded by the retained
tube/projections. Therefore the innovation sequence required to maintain
IE-1 is itself uniformly bounded. There is no contradiction from innovation
blow-up.

The literal mean update is a bounded-input finite-dimensional recursion:
prediction has stable AW OU factor, BA OU factor when enabled, bounded
attitude transport, and the S pseudo-update supplies feedback to the
integrated chain; accepted accelerometer/magnetic gains are finite on the
retained covariance class. A uniformly bounded innovation sequence therefore
does not force the estimator mean to leave the retained set by any currently
proved inequality. In particular there is no cumulative innovation-energy
budget in MARINE MOTION/IMU BIAS/MAGNETIC SERVICE.

Thus the q->0 limiting equations are dynamically **compatible with bounded
base inputs at the level of all proved bounds**. This rules out an exclusion
proof based only on boundedness, OU decay, S feedback, or physical jerk.

However it still does not construct an infinite shipping execution: the
bounded innovations must be generated by one globally admissible physical
trajectory while the estimator recursion remains exactly on the compatibility
manifold and MAGNETIC SERVICE remains strict for all time. The local IFT
construction proves finite-horizon continuation near a suitable strict-margin
base, but no global continuation theorem prevents loss of transverse gain
rank or service margin over infinitely many words.

Therefore the infinite-equality problem reduces to **global continuation of
the local compatibility manifold**. A decisive theorem must prove one of:

(A) continuation: the compatible manifold contains a forward-complete
strict-margin A21 trajectory; then an infinite equality execution exists and
O2 block contraction is false;

(B) escape: every compatible forward trajectory reaches in finite time a
boundary where compatibility, transverse control rank, retained-state
bounds, or MAGNETIC SERVICE fails; compactness then yields finite m,delta and
O2 block contraction.

No current invariant proves A or B. The existing equations show only that
there is no elementary boundedness contradiction as q->0. Hence the infinite
execution question remains open, but its precise mathematical form is now a
forward-completeness/escape problem for the closed-loop compatibility
manifold, not an observability or covariance problem.
## Viability criterion for the exact compatibility manifold

Let h(m,q)=0 denote the two transverse literal accelerometer compatibility
constraints at an applied correction, with m the base estimator/physical
state and q the BA/tilt ratio. Between corrections q contracts by the known
BA factor and m follows the shipping prediction/S/magnetic maps. At an
accelerometer correction the control is the physical/base innovation u.

The discrete viability condition is

`h(F_k(m,u), Phi_b q)=0`.                                   (VIAB-1)

Linearizing in the two transverse control components gives the tangency
matrix
`G_perp = d_u_perp h(F_k(m,u),Phi_b q)`.                    (VIAB-2)
If `sigma_min(G_perp)>=g_c>0` and the required control remains inside strict
physical/gate margins, the implicit-function theorem gives a unique local
compatible continuation. Thus loss of compatibility can occur only by
reaching one of four boundaries: `det G_perp=0`, physical/retained-state
boundary, event/gate boundary, or magnetic-service boundary.

In the q->0 limit the compatibility target h=0 becomes nominal field-axis
collinearity. The required correction per step remains bounded (indeed the BA
target change vanishes with q); no term in VIAB-1 diverges. Therefore there
is no intrinsic vector-field singularity forcing escape as q->0.

Conversely, forward completeness would follow from a compact positively
invariant subset of the compatibility manifold on which all four margins are
strictly positive. The current proof has no lower bound on `sigma_min(G_perp)`
for reachable A21 states and no invariant lower margin for MAGNETIC SERVICE
or gates along the constrained continuation. Hence such a subset is not
proved to exist.

An escape theorem would require a scalar function V_c on the compatibility
manifold whose increment has a strict sign until a boundary is hit. BA norm
`|q|` cannot serve: it decreases toward the interior zero-BA manifold rather
than a forbidden boundary. The retained estimator energies and S-chain
storage are actively stabilized and likewise have no proved monotone drift
toward a boundary. Physical bounded-potential motion is not a state function
of the nominal compatibility manifold.

Therefore neither forward completeness nor finite escape follows from the
current assumptions. More strongly, the local equations show q=0 is a
regular candidate limiting manifold rather than a forbidden boundary, so an
escape proof cannot be based on BA decay alone.

**Resolution status.** The global continuation/escape problem is undecidable
from the currently proved invariants: the necessary viability rank/margin
conditions are not controlled source-uniformly, and no monotone escape
functional exists in the theorem set. Proving either branch requires a new
shipping invariant (uniform compatibility-control rank/service margin) or an
explicit forward-complete construction. This is the exact O2 theorem gap.
## Full OU/S/accelerometer adaptation zero dynamics

Use the literal shipping order: OU/LIN prediction, pending AW covariance
inflation, scheduled S=0 pseudo-update inside prediction, then accelerometer
update; AW covariance sync/inflation is mean-neutral. Let m collect the
transverse nominal mean blocks `(a_w,v,p,S,b_a,theta,...)` immediately after
an accelerometer correction.

For one prediction interval write the exact mean prediction as
`m^- = F_k m^+`, with the shipping OU/LIN/BA/attitude coefficients. If an S
update is due,

`m^S = (I-K_S,k H_S) m^-`,                                  (ZD-1)

because `r_S=-H_S m^-`. This is endogenous feedback, not an external input.
At the following accelerometer correction,

`m^(+)next = R_k[ m^S + K_a,k r_a,k ]`,                     (ZD-2)

where R_k denotes the literal quaternion/error reset map on the mean
coordinates (identity on unaffected Euclidean blocks). Covariance sync does
not enter ZD-1--ZD-2 at mean level, though it changes future gains.

Define the transverse compatibility output at that accelerometer epoch
`e_c,k=P_b,k fhat_k=P_b,k C_f,k m^S` in the q->0 limit. Exact persistent
compatibility imposes

`e_c,k=0`.                                                   (ZD-3)

The physical measurement identity gives
`r_a,k=f_phys,k-fhat_k-bhat_a,k`.                           (ZD-4)
Substituting ZD-4 into ZD-2 and then the next prediction/S update yields the
closed forced zero-dynamics map

`m_(k+1)^S = A_cl,k m_k^S + B_phys,k f_phys,k`,             (ZD-5)

with
`A_cl,k=(I-K_S H_S) F_k R_k [I-K_a(C_f+C_ba)]`
(with the literal ordering/frame blocks), and
`B_phys,k=(I-K_S H_S) F_k R_k K_a`.                        (ZD-6)

The constraint `P_b C_f m_k^S=0` is imposed at every applied accelerometer
epoch. Thus S feedback is fully inside A_cl: its memory cannot be chosen
independently.

### Does S feedback exclude an infinite zero trajectory?

No source-uniform exclusion follows from the current coefficients. ZD-5 is a
bounded linear/time-varying forced recursion on the retained tube. The
physical force is an admissible bounded forcing, and S feedback is stabilizing
rather than an accumulating conserved quantity. Setting e_c=0 imposes two
linear constraints per epoch on the forced recursion, while f_phys has three
components plus the independently evolving admissible attitude history.

To prove impossibility one would need a left-annihilator of the constrained
input map: a nonzero row L_k such that `L_k B_phys,k=0` but the zero-output
condition forces `L_k A_cl,k m` to a nonzero value under recurring MARINE
MOTION. Equivalently, the constrained Rosenbrock matrix of the literal
adaptation cycle would need a source-uniform rank condition excluding an
invariant zero.

The current proof has no such rank theorem. In fact the previously proved
local rank-two transverse AW correction at diagonal covariance indicates the
opposite locally: where the relevant gain block has rank two, the physical
accelerometer forcing has enough transverse authority to enforce ZD-3 while
the S feedback is simply part of the state transition. Endogenous S memory
changes the required forcing but does not by itself overdetermine it.

Therefore incorporating the S pseudo-measurement correctly does **not yet**
exclude persistent nominal field-axis compatibility. It sharpens the missing
lemma to an invariant-zero/rank statement for `(A_cl,B_phys,C_c)` over the
actual gain/covariance history. A proof that this constrained system has no
bounded invariant zero under recurring MARINE MOTION would close O2; a
forward-complete bounded zero-dynamics solution would refute block
contraction.
## Literal lifted invariant-zero formulation (terminology corrected)

`applyIntegralZeroPseudoMeas()` is an ordinary Kalman pseudo-measurement
update with measurement `S=0` and innovation `r_S=-S`. There is no separate
controller or extra adaptation loop. References below to "S feedback" mean
only the state-dependent effect of this standard measurement update through
the Kalman cross-covariance gain; use "S pseudo-update" in theorem text.

Let k_j be successive scheduled S pseudo-update epochs. Lift the literal
shipping recursion from immediately after S update j to immediately after S
update j+1. Between them compose, in actual chronological order, every OU/LIN
prediction, accelerometer update, magnetic update, quaternion reset, BA
projection and covariance-sync operation. Covariance sync changes future
gains but has no direct mean increment.

For a fixed actual base history H_j over this interval, the first-order mean
map has the affine form

`m_(j+1)=A_j(H_j)m_j + B_j(H_j) u_j + d_j(H_j)`,             (LZ-1)

where u_j stacks the physical accelerometer/magnetic measurement forcing over
the interval. The gains in A_j,B_j are the literal gains generated by the
same covariance history; they are not independent design variables.

Stack every applied accelerometer compatibility output in the interval:

`y_j=C_j(H_j)m_j + D_j(H_j)u_j + e_j(H_j)`.                 (LZ-2)

In the q->0 equality limit, persistent compatibility requires `y_j=0` for
all j. A bounded invariant-zero execution is therefore a same-history
sequence satisfying LZ-1--LZ-2, all shipping gates/projections, MARINE MOTION
and MAGNETIC SERVICE, with bounded m_j.

Eliminate u_j only if the literal constrained input block D_j has the
required rank. The finite-interval Rosenbrock matrix is

`R_j(z)=[[zI-A_j,-B_j],[C_j,D_j]]`,                         (LZ-3)

but because coefficients are time varying, absence of a zero of one frozen
R_j is neither necessary nor sufficient. The needed theorem is a lifted
uniform left-invertibility/zero-dynamics statement over admissible sequences.

Existing results do not prove it. In particular, the S pseudo-update adds
rows/couplings to A_j through `K_S(-S)`, but its residual is endogenous and
does not reduce the dimension of physical forcing by itself. Where D_j has
full transverse row rank, y_j=0 can be solved locally for transverse physical
forcing; the remaining question is whether the resulting constrained lifted
map admits a bounded forward-complete same-history sequence.

Thus the exact invariant-zero blocker survives correct treatment of the S
pseudo-update. What has improved is the formulation: it is now the lifted
time-varying shipping map LZ-1--LZ-2, with gains/covariance chronology
included exactly. No claim that a separate S feedback controller exists is
made.
## Long-horizon invariant-zero test against bounded physical motion

Consider a bounded forward-complete zero-output execution of LZ-1--LZ-2 in
the q->0 regime. At applied accelerometer epochs the nominal specific force
is field-axis compatible. The physical forcing is
`f_phys=fhat+bhat_a+r_acc`.                                 (IZ-1)

Sum the physical translational dynamics over N lifted intervals. MARINE
MOTION's bounded velocity implies the time average of physical acceleration
tends to zero along long horizons; bounded displacement/potential further
excludes a persistent nonzero DC component of velocity/position. Therefore
an invariant-zero construction would be impossible if compatibility required
a source-uniform nonzero DC transverse component of physical acceleration.

But the literal zero-output constraint does not require such a DC component.
The field axis and attitude may rotate under admitted MOVING motion, and the
required transverse accelerometer innovation can be oscillatory with zero
mean. The S pseudo-update acts on the estimator integral state and can also
produce bounded periodic corrections. Thus all long-horizon sums can cancel
without violating bounded velocity/displacement/potential.

Indeed the zero-output equations are affine pointwise constraints on the
physical measurement forcing, not a sign-definite work/energy identity.
Neither the Kalman gains nor the S pseudo-update make the required physical
forcing have a fixed sign or nonzero mean. Consequently summation/telescoping
of LZ-1 does not yield a contradiction with MARINE MOTION.

### Theorem consequence

Under the current assumptions, a bounded forward-complete invariant zero
cannot be **excluded** by any existing long-horizon physical boundedness
condition. The remaining forcing can be zero-mean and periodic/quasiperiodic.
However existence is still not proved because one must solve the coupled
time-varying gain/covariance/mean recursion on the zero-output manifold.

A constructive existence theorem would require a periodic (or compact
recurrent) fixed point of the lifted constrained Poincare map. The natural
candidate is a periodic MARINE MOTION attitude/force history with strict
magnetic service and zero-mean translational acceleration, together with a
periodic covariance/tuner orbit. No analytical fixed-point theorem for that
shipping map is currently established.

Therefore the invariant-zero question remains undecided, but one branch is
now narrowed further: **finite escape cannot be proved from bounded physical
velocity/displacement/potential or the S pseudo-update alone**, because the
required zero-output forcing need not contain a forbidden DC component.
## Reduced zero dynamics after eliminating transverse physical forcing

On a fixed lifted interval write the compatibility output as
`0=y=Cm+D_perp u_perp+D_parallel u_parallel+e`.              (RZ-1)
Whenever the literal transverse input block D_perp is invertible, eliminate

`u_perp=-D_perp^-1(Cm+D_parallel u_parallel+e)`.             (RZ-2)

Substitution into the lifted mean map gives the exact constrained dynamics

`m_+=A_z m+B_z u_parallel+d_z`,                             (RZ-3)

`A_z=A-B_perp D_perp^-1 C`,
`B_z=B_parallel-B_perp D_perp^-1 D_parallel`.               (RZ-4)

All S pseudo-update effects are already inside A,B,C,D through the literal
chronology and gains. Thus persistence reduces locally to boundedness and
admissibility of this one-input time-varying zero dynamics.

Project RZ-3 onto the LIN chain `(v,p,S)` and AW state. The open-loop
prediction has the polynomial-integrator chain driven by exponentially
stable AW; each scheduled S pseudo-update is a Kalman measurement correction
with innovation `-S`. Consequently the lifted S-to-S map has no conserved
integrator mode forced by the zero-output constraint: its homogeneous
coefficients are finite, while the remaining longitudinal input can alter the
chain through the accelerometer gain/cross covariance.

Critically, the compatibility output removes only the two transverse force
components. It imposes no condition on the longitudinal physical-force
component. That remaining scalar input is sufficient to cancel a scalar
secular S/velocity condition if one appears; there is no second independent
longitudinal closure equation in the persistent-pair problem.

Therefore the literal S pseudo-update does not produce an algebraic
overdetermination of the reduced zero dynamics. A finite-escape theorem would
require showing that RZ-3 has an unstable/unbounded mode that is both
uncontrollable from B_z and unavoidable under MARINE MOTION. No such mode is
present in the current zero-action classification: four-S injectivity removes
homogeneous LIN/AW zero-action roots, while the base forced recursion may
remain bounded.

This means the invariant-zero analysis does not currently yield an exclusion
theorem. Locally, where D_perp stays nonsingular and margins are strict, the
zero manifold is viable and the remaining longitudinal input can maintain
bounded integral states. Global existence still depends on preventing loss of
D_perp rank/margins, but there is no intrinsic S-chain escape mechanism.
## Global controlled-invariance test for the reduced zero dynamics

After transverse elimination the local compatibility dynamics are
`m_+=A_z(m,H)m+B_z(m,H)u_parallel+d_z(m,H)`, with the actual covariance/
gain/history variables H carried in the augmented state. Let K be a compact
strict-margin subset on which D_perp is nonsingular and all shipping branches
are fixed.

A forward-complete invariant zero exists if K contains a nonempty controlled
viability kernel: for every state in that kernel there is an admissible
longitudinal input with the next augmented state again in the kernel. Finite
escape follows if every candidate compact K has a boundary point whose
outward normal n satisfies `n'B_z=0` and `n'(A_zm+d_z-m)>0` uniformly (or the
corresponding discrete tangent-cone condition).

The literal shipping structure supplies neither condition source-uniformly.
`B_z` is one-dimensional and cannot control arbitrary normals of the full
augmented mean/covariance/tuner state, so a generic controlled-invariance
theorem is unavailable. Conversely, covariance/tuner evolution is autonomous
given the event/mean history and is bounded by existing guards; no proved
boundary normal has an unavoidable outward drift. S, AW and BA mean dynamics
are dissipative/bounded rather than sign-definitely escaping.

Most importantly, D_perp depends on covariance cross blocks and the next
accelerometer Jacobian. Previous analysis proved neither a positive reachable
lower bound on sigma_min(D_perp) nor inevitable approach to zero. MAGNETIC
SERVICE likewise has a closed lower floor for admissible words but no theorem
that its strict margin monotonically decreases along constrained zero
dynamics.

Therefore the global question cannot be decided by controlled-invariance
geometry from the current invariants: there is no certified compact viable
set and no certified escape boundary. The exact undecided quantities are now
`sigma_min(D_perp)` and the strict gate/service margins along the constrained
trajectory.

This is a genuine assumption/theorem gap, not a missing algebraic manipulation.
To prove O2 under current assumptions one must derive an additional invariant
from shipping showing finite-time loss of one of those quantities. To refute
O2 one must construct a forward-complete constrained trajectory with their
positive infima. Neither follows from the presently proved bounds.
## Compatibility at the exact MAGNETIC SERVICE boundary

MAGNETIC SERVICE is the closed condition
`lambda_min(G_M)>=mu_M`, with `mu_M>0`, formed from actually applied
transported magnetic rows on the service sector. At the boundary
`lambda_min(G_M)=mu_M`, the service Gram remains strictly positive definite
on that sector. No magnetic information direction is lost merely because the
surplus above mu_M vanishes.

The surviving physical compatibility line is not a null direction of this
service Gram in the same coordinate space. Magnetic rows first restrict the
attitude/gyro root to the transported field-axis compatibility class; the
one-dimensional tilt/BA line survives only after combining that magnetic
restriction with accelerometer/BA compatibility. Thus service equality does
not force the physical compatibility generator to disappear.

All literal magnetic Joseph updates and chronological transports are
continuous at `lambda_min(G_M)=mu_M`; their innovation covariances retain the
positive Rmag floor. Therefore a same-history compatible sequence can converge
to a service-boundary word while retaining well-defined magnetic corrections
and the one-dimensional compatibility relation.

Consequently there is no analytical escape mechanism at the exact service
boundary. Compatibility can persist **at the level of the theorem
constraints** when service equals mu_M; the assumption intentionally admits
that boundary. A global continuation proof does not need a strict service
surplus if it can remain on/inside the closed service set.

This removes magnetic-service surplus as an intrinsic O2 blocker. What remains
is whether the constrained shipping dynamics can keep the actual service
inequality >=mu_M while maintaining compatibility; equality itself is not a
singularity. No current invariant forces crossing below mu_M.

Hence the only local viability quantity capable of destroying compatibility
before violating an explicit physical assumption is the transverse
compatibility-control rank `sigma_min(D_perp)`. The O2 global problem reduces
to whether D_perp can remain nonsingular along an admissible constrained
trajectory (with service allowed to sit on its closed boundary).
## Can the shipping recursion maintain transverse compatibility authority?

Let D_perp,k be the derivative of the next compatibility output with respect
to the two transverse components of the physical accelerometer measurement at
an accepted correction. Because the measurement enters the mean only through
`K_a r_a`, D_perp factors as

`D_perp,k = C_(c,k+1) T_(k+1<-k) K_a,k E_perp,k`,           (DA-1)

including the literal reset/prediction/S-pseudo-update transport T and the
next compatibility-output derivative C_c. Since
`K_a=P C_a' S_a^-1`, loss of authority can arise from the projected gain
numerator or from its terminal/output transport.

Prediction injects fresh AW process covariance
`Q_aa=(1-phi_a^2) Sigma_aw_stat >0` into the AW marginal without an equal
fresh AW-attitude/BA cancellation term. Thus immediately after every positive
prediction interval the accelerometer gain numerator contains a favorable
fresh term `Q_aa R_wb'`. However its magnitude is O(h/tau), while inherited
cross-covariances remain O(1). No Loewner/sign invariant prevents exact
cancellation of the **projected** numerator after adding the inherited terms.

Repeated accepted accelerometer updates do not force rank loss: for a fixed
accelerometer Jacobian their own covariance update right-multiplies `P C_a'`
by an invertible 3x3 factor, preserving its row-block rank. Prediction adds
fresh AW variance, and reset is invertible. Therefore there is likewise no
structural mechanism that necessarily drives D_perp toward singularity.

S and magnetic pseudo/measurement updates can change the relevant
cross-covariances, and the next C_c/J_att geometry changes with the nominal
mean. Their covariance maps are continuous bounded Joseph maps with positive
measurement-noise floors. They have no monotone determinant law for
D_perp. Consequently the shipping recursion has neither a proved repelling
barrier from `det D_perp=0` nor an attracting mechanism toward it.

### Structural conclusion

The actual covariance recursion is compatible with maintaining D_perp
nonsingular indefinitely: none of its operations forces rank loss, and
prediction recurrently injects favorable AW variance. But this is not an
existence theorem for a particular forward-complete compatible execution,
because inherited cross-covariances and changing output geometry could still
drive a chosen trajectory to cancellation.

Conversely, a theorem that every compatible trajectory must reach
`det D_perp=0` is impossible from the presently known operation-wise
identities: accelerometer self-updates preserve rank and prediction can move
away from cancellation. Such a theorem would require a new global invariant
coupling S/magnetic cross-covariance evolution to compatibility; none is
present in shipping/proof assumptions.

Therefore rank loss of D_perp cannot serve as a currently proved mandatory
escape mechanism. At the structural level indefinite nonsingularity is
allowed, but actual existence remains a global reachability question.
## Near-constant world field closes the q->0 geometry to a DC AW requirement

MAGNETIC SERVICE now includes a physical world-frame geomagnetic reference
`b_M^W`, field variation `||b^W-b_M^W||<=eps_B`, and a uniform gravity/field
non-collinearity margin `g_Bperp>0`.

In the q->0 exact compatibility limit, the accelerometer attitude row
`J_att=-[f_cog^B]x` annihilates the transported magnetic-axis tilt. Therefore
`f_cog^B` is parallel to the body magnetic direction. Rotating both vectors
to world coordinates preserves collinearity. Since

`f_cog^B=R_wb (a_w^W-g^W)`,

the exact constant-field case implies

`P_(b_M^W)^perp (a_w^W-g^W)=0`, hence
`P_(b_M^W)^perp a_w^W=P_(b_M^W)^perp g^W`.                 (GF-1)

With near-constant field and bounded lever/reference defects,

`||P_(b_M)^perp a_w^W-P_(b_M)^perp g^W|| <= C_B eps_B+C_L`, (GF-2)

where C_L contains only the explicitly bounded lever/linearization terms.
Thus if `C_B eps_B+C_L<g_Bperp`, persistent q->0 compatibility requires a
nonzero world-fixed transverse component of nominal AW bounded below by

`a_DC := g_min cos(80 deg)-C_B eps_B-C_L >0`.                  (GF-3)

Numerically `cos(80 deg)=0.1736481777`, so before the explicit field-variation/lever deductions the compatibility mode requires at least `0.1736481777 g_min` of world-fixed nominal AW transverse to the geomagnetic direction. For gravity near 9.81 m/s^2 this is about 1.70 m/s^2; this is a gravity-scale requirement, not a noise-floor trajectory.

This is qualitatively different from an arbitrary oscillatory retuning: the
required component has a fixed world direction inherited from gravity and the
near-constant Earth field.

### Interaction with the LIN/S chain

Between measurement corrections AW is OU-decayed, but accelerometer and S
pseudo-updates can change its mean. If GF-1/GF-3 holds at every dense applied
accelerometer epoch indefinitely, the nominal AW samples contain a persistent
world-fixed DC component. The deterministic LIN chain integrates AW into v,p,S.

To turn this into an exclusion theorem one must account for measurement
corrections to v,p,S as well as AW: the S pseudo-update and accelerometer
cross-gains can remove integrated drift. Therefore the bare identity
`v_dot=a_w` between updates is insufficient to claim unbounded v/p/S.

The new field premise nevertheless removes the previous zero-mean forcing
escape: any persistent compatibility construction must now use the Kalman
measurement corrections themselves to cancel the fixed DC LIN injection.
The remaining exact lemma is finite-dimensional and filter-specific:

`persistent GF-3 + bounded (v,p,S)`
` => measurement-correction action has a nonzero DC component of at least
   c_DC(a_DC)>0` on every sufficiently long window.          (GF-4 target)

If the literal correction maps cannot supply that DC cancellation while
remaining on the compatibility manifold, infinite equality is excluded. If
they can, that supplies the explicit persistent-zero construction mechanism.

Thus near-constant geomagnetism materially sharpens O2 but does not by itself
complete it; the remaining blocker is now the DC balance of the **literal
measurement-corrected LIN/S mean recursion**, not arbitrary physical-to-
nominal tracking.
## DC balance of the literal measurement-corrected AW/LIN chain

Let `e_B` be any fixed unit world direction in the gravity component
transverse to the near-constant geomagnetic field. Persistent q->0
compatibility gives, up to the explicit eps_B/lever defect,
`e_B' a_w,k >= a_DC>0` at every dense accelerometer epoch.

Between corrections the OU mean prediction is
`a_w^- = phi_a a_w^+`. Therefore each positive prediction interval removes
`(1-phi_a)e_B'a_w^+` from this fixed component. Summing the exact AW mean
recursion over N steps gives

`a_w,N-a_w,0
 = -sum_k (1-phi_a,k)a_w,k^+
   +sum_(S updates l) K_awS,l (-S_l^-)
   +sum_(acc updates k) K_awa,k r_a,k
   +other literal mean-reset terms`.                        (DCB-1)

Quaternion reset does not directly change the world AW state; covariance sync
is mean-neutral. For bounded a_w and persistent `e_B'a_w>=a_DC`, division by
elapsed time yields the necessary average balance

`liminf (1/T) e_B' [sum_S -K_awS S + sum_a K_awa r_a]
` >= lambda_a a_DC - defects >0`,                           (DCB-2)

where `lambda_a` is the source-uniform OU decay rate implied by the tau/dt
bounds.

Thus an infinite equality execution requires a nonzero DC **Kalman correction
supply** into AW. This is now rigorous at the mean-recursion level.

Can that supply occur while physical acceleration has zero long-time mean?
Yes in principle: `r_a=f_phys-fhat-bhat_a`. Under compatibility, fhat contains
the gravity-canceling nominal AW component, so even zero-mean physical
acceleration produces a residual with a nonzero DC component relative to the
misaligned nominal prediction. The accelerometer update can therefore supply
the positive average in DCB-2. The S pseudo-update may add or subtract from
the same balance. No current assumption bounds the long-time mean of
`K_awa r_a` to zero.

Consequently the new geomagnetic floor proves that persistence requires
continuous nonzero correction work, but it does **not** make that work
unavailable: the physical accelerometer residual generated by the erroneous
DC nominal AW supplies exactly such a correction channel.

To obtain a contradiction one would need a gain identity showing that the
same DC residual necessarily drives the nominal AW **away** from, rather than
back toward, the compatibility value after accounting for S updates. The
Kalman sign is the opposite: for the direct AW measurement block
`J_aw=R_wb`, the correction tends to reduce the physical-vs-nominal force
residual. Thus the DC balance is compatible with a steady biased nominal AW
only if cross-coupled S/attitude/BA corrections continually recreate the
bias. Whether that closed balance has a fixed point is the remaining exact
question.
## Sign of the fixed-direction DC adaptation loop

Project the mean dynamics onto the fixed world direction e_B from GF-3. Let
`a=e_B'a_w`. Ignore only terms already carried as explicit bounded defects.

**OU prediction.** `a^- = phi_a a^+`, with `0<phi_a<1`: strict attenuation.

**S pseudo-update.** Its AW increment is
`Delta a_S = - e_B' K_awS S`.                               (LG-1)
The sign of this term is not fixed by PSD covariance. `K_awS=P_aw,S
(P_SS+R_S)^-1`; the cross covariance P_aw,S may have either sign along e_B.
Thus the S pseudo-update is not intrinsically negative or positive feedback
on a DC AW component.

**Accelerometer update.** In measurement coordinates the direct AW-only
conditional factor is
`I-K_aw R_wb`, with `K_aw=P_aw,C S_a^-1`. If all cross covariances with
attitude/BA/gyro were absent, `P_aw,C=P_aw,aw R_wb'` and the scalar gain
`k_a=e_B' K_aw R_wb e_B` lies in (0,1), so a physical zero-mean acceleration
would contract a nominal DC AW bias by factor `1-k_a`.

With the literal filter, however,
`P_aw,C=P_aw,theta J_att' + P_aw,aw R_wb' + P_aw,ba + ...`. (LG-2)
The cross terms can change both magnitude and sign of the projected scalar
`k_a=e_B'K_aw R_wb e_B`. PSD of P and S_a>0 do not imply `0<k_a<1` for this
sub-block. The full Kalman covariance contracts measurement uncertainty, not
each individual state component monotonically.

Therefore the one-cycle homogeneous DC coefficient has schematic exact form

`L_DC = (accelerometer/reset transport)
`       * (I-K_S H_S) * phi_a`,                             (LG-3)

projected from e_B-AW back to e_B-AW, plus cross-state return paths through
`v,p,S,theta,ba`. Those return paths are part of the same matrix product and
can have either sign.

### Result

No sign theorem `0<=L_DC<1` follows from Kalman/Joseph structure alone.
Likewise equality `L_DC=1` is not structurally guaranteed. The complete DC
loop gain is history/cross-covariance dependent.

The exact condition for a nonzero DC compatibility fixed point is
`(I-L_DC) a_DC = d_DC`,                                    (LG-4)
where d_DC is the projected constant forcing generated by gravity/reference
and physical measurement terms after imposing compatibility. If d_DC=0,
a nonzero homogeneous fixed point requires eigenvalue one; with nonzero
forcing a bounded biased fixed point may exist even when |L_DC|<1.

This corrects the proposed dichotomy: proving total loop gain <1 does **not**
by itself exclude a nonzero required a_DC; a stable affine loop can support a
nonzero fixed point. Exclusion requires showing the unique affine fixed point
is inconsistent with GF-3/physical measurements, or that the required forcing
cannot be supplied.

Thus the next exact calculation is the affine fixed-point equation of the
complete lifted DC subsystem, not merely the sign of its homogeneous gain.
## Affine DC fixed-point equation for the lifted adaptation state

Project onto the fixed world direction e_B and collect every mean coordinate
that can return into AW under one lifted scheduler pattern:
`x_D=(a_w,v,p,S,b_a,theta,...)_eB`. For a fixed repeated same-history
coefficient pattern H and admissible constant/periodic physical forcing, the
literal lifted map is affine

`x_D^+ = A_D(H) x_D + d_D(H,u_phys)`.                       (AF-1)

`A_D` is the exact product of prediction, scheduled S pseudo-update,
accelerometer/magnetic Kalman updates and resets; covariance sync affects A_D
through the gains but has no additive mean term.

A periodic DC fixed point satisfies

`(I-A_D)x_* = d_D`.                                         (AF-2)

If I-A_D is invertible,
`x_*=(I-A_D)^-1 d_D`; otherwise solvability requires the Fredholm condition
`l'd_D=0` for every left unit eigenvector `l'A_D=l'`, with neutral components
then fixed by compatibility/boundedness.

The geomagnetic q->0 compatibility condition adds

`C_B x_* = g_Bperp + delta_B`,                              (AF-3)

where C_B selects the world AW component and
`|delta_B|<=C_Bfield eps_B+C_L`; in the exact constant/no-lever case
`C_B x_*=g_Bperp`.

Combining AF-2--AF-3 gives the exact algebraic solvability condition

`C_B (I-A_D)^-1 d_D = g_Bperp+delta_B`                      (AF-4)

when invertible, together with all remaining compatibility rows. This is the
affine invariant-zero fixed-point equation.

Now expose d_D. The physical accelerometer enters only through innovations
`r=f_phys-fhat-bhat_a`; after substitution, d_D contains the physical
specific-force/gravity/reference terms multiplied by the literal Kalman gains.
Therefore d_D is not zero even for zero-mean translational acceleration:
gravity and the imposed nominal compatibility offset appear in the affine
forcing.

Consequently no source-independent contradiction follows from AF-4. Depending
on the actual gains/cross covariances, the affine map can in principle have a
fixed point with the required `C_B x_*=g_Bperp`. Conversely, no theorem says
it must.

For nonperiodic admissible marine motion the corresponding exact statement is
the bounded particular-solution equation
`x_k=Phi_A(k,0)x_0+sum_i Phi_A(k,i+1)d_i`,                  (AF-5)
with `C_B x_k` constrained by GF-2 at every accelerometer epoch. Periodic
AF-4 is therefore a sufficient construction mechanism, not a necessary form
of every persistent zero.

Thus the affine calculation identifies the decisive quantity:
`C_B(I-A_D)^-1 d_D` for a candidate recurrent shipping pattern. Proving it
uniformly separated from the geomagnetic target would exclude periodic
persistent zeros; exhibiting equality with strict margins would construct
one. Current symbolic shipping bounds do not determine its value or sign.
## Marine-motion window decomposition for the geomagnetic O2 argument

The constant/repeated-forcing affine fixed point is only a possible
counterexample construction device, not the controlling exclusion proof.
MOVING retains its full time-varying physical history.

On every excitation window W=[t,t+T_E], project onto the fixed world direction
e_B supplied by the near-constant geomagnetic premise and decompose nominal
and physical acceleration into a window mean plus a zero-mean fluctuation.
Exact q->0 compatibility gives a positive lower bound a_DC on the nominal AW
window mean in this direction, while the fluctuation retains the actual marine
time variation.

Physical acceleration is not set to zero or constant. Its integral is the
physical velocity increment, and concatenated windows retain all MARINE
MOTION velocity/displacement/potential, acceleration/jerk, and recurring
attitude-span conditions.

Summing the literal measurement-corrected AW recursion over W gives two
coupled requirements: (1) the correction sequence must replace the positive
mean OU loss associated with a_DC; (2) the same corrections must accommodate
the zero-mean time-varying marine forcing while the LIN/S states remain
bounded. Any exclusion theorem must hold uniformly over these admissible
fluctuation histories. Absence of a constant affine fixed point would not
suffice.

### Exact varying-history compatibility identity

At every applied accelerometer epoch k in the q->0 limit,

`P_Bperp [a_w,k^W-g^W]=delta_B,k`,                         (MM-1)

with `||delta_B,k||` bounded by the explicit near-field/lever defects.
Thus the transverse nominal AW is pinned near the fixed world vector
`P_Bperp g`; it is not free to follow the varying physical acceleration.

Using the exact innovation identity and rotating to world coordinates,

`P_Bperp R_bw r_a,k
 =P_Bperp a_phys,k^W-P_Bperp a_w,k^W
  -P_Bperp R_bw bhat_a,k + lever/calibration terms`.        (MM-2)

Substituting MM-1,

`P_Bperp R_bw r_a,k
 =P_Bperp a_phys,k^W-P_Bperp g^W
  -P_Bperp R_bw bhat_a,k + bounded defects`.                (MM-3)

Hence recurring MARINE MOTION variation appears explicitly in the
accelerometer residual around a fixed gravity-sized offset. The residual is
not an arbitrary control.

For persistence, the complete S/accelerometer correction sequence must map
the varying residuals MM-3 back to states satisfying MM-1 at every subsequent
accelerometer epoch. This is the correct filter-specific bridge. The remaining
lemma is to show that the literal lifted correction map cannot annihilate all
allowed zero-mean physical variations while also supplying the fixed positive
OU-loss replacement and keeping LIN/S bounded. This must use the MARINE
MOTION excitation span; a DC-only fixed-point argument is insufficient.

## Rotating accelerometer-innovation action from gravity span and geomagnetism

The actual MOVING hypothesis is the stronger/specific gravity-direction span
condition `Delta_g(W)>=theta_E` on every complete moving excitation window,
not an unspecified attitude norm span.

Let `b` be the near-constant unit world geomagnetic direction and decompose
`g=g_parallel+ d`, with `d=P_bperp g`. The inclination restriction gives
`|d|>=g_min cos(80 deg)`.

In the q->0 compatibility limit the nominal CoG specific-force vector is
parallel to b, so after subtracting the physical measurement model the
body-frame accelerometer innovation contains the rotating term `R_wb(t)d`
(up to physical translational acceleration, BA and the explicit bounded
defects). Because rotations preserve norm, this term has fixed magnitude at
least the geomagnetic floor.

Gravity span alone does not imply a positive span of `R_wb d`: rotations
about d are a geometric exception. This loophole must not be hidden. Combine
gravity span with the magnetic service geometry. Since g and b are fixed
non-collinear world vectors, the map R -> (Rg,Rb) is injective on SO(3).
On the compact set with angle(g,b) bounded away from zero by the inclination
cutoff, there is a modulus c_gb>0 such that any orientation pair with gravity
span at least theta_E produces a positive joint span of the two body vectors.
If `R d` had zero span while `R g` had positive span, the variation must
occur through the b component; the applied magnetic rows then carry that
variation. Thus the **joint accelerometer+magnetic signed record**, rather
than accelerometer innovation alone, has a source-uniform rotating action.

Formally define on W
`A_rot(W)=sum_acc ||P_record R_wb,k d||^2 + sum_mag ||Delta(R_wb,k b)||^2`
with the literal whitening/transport used by the complete word. Compactness
of SO(3), `|d|>=g_min cos80`, `Delta_g>=theta_E`, and MAGNETIC SERVICE give
a positive existence modulus
`A_rot(W)>=a_rot(theta_E,mu_M,I_max,eps_B)>0`.

This is an existence-level geometric modulus; an explicit closed formula for
a_rot is not yet derived. It avoids the false claim that gravity span alone
forces R d to rotate.

The next adaptation step must retain the **signed joint record**. Its
accelerometer part enters the literal Kalman mean update, while its magnetic
part enters the attitude correction/reset and thereby changes the subsequent
accelerometer Jacobian. The target is to show that an infinite q->0
compatibility execution would have to annihilate this recurring positive joint
action through the same-history adaptation map, contradicting the complete-word
nullity classification unless the physical tilt/BA kernel is nonzero; as
q->0 that kernel has zero BA component and the geomagnetic floor forbids the
remaining pure tilt compatibility.

## Zero-BA endpoint lemma and neighborhood consequence

Consider a complete admissible MOVING word and a zero-action slow-root vector
with q=0. The existing four-S and magnetic/gyro arguments first remove LIN/AW
root and gyro-bias quotient components. The remaining candidate is pure
attitude tilt theta.

At every applied magnetic row, zero loss gives
`[R_wb b^W]x theta_k=0`, so the transported nonzero theta_k must be parallel
to the body magnetic direction. At every applied accelerometer row, zero loss
with q=0 gives
`[R_wb(a_w^W-g^W)]x theta_k=0`, so the same transported tilt must be parallel
to the nominal CoG specific-force direction. Hence a nonzero pure-tilt null
trajectory requires
`a_w^W-g^W parallel b^W` at every relevant accelerometer epoch.          (ZB-1)

Under the near-constant nonvertical geomagnetic premise this implies
`P_bperp a_w^W=P_bperp g^W+delta_B`, with the positive gravity-transverse
floor from GF-3.

However, ZB-1 by itself is a condition on the **nominal** AW mean. The
complete-word homogeneous nullspace equations do not constrain the base AW
mean through MARINE MOTION; physical gravity span constrains the true
orientation/history, while accelerometer innovations may make the nominal AW
satisfy ZB-1. Therefore the stated assumptions do **not** by themselves imply
that a pure-tilt homogeneous null is impossible.

This exposes an error in the proposed endpoint shortcut: the joint rotating
base innovation action is not part of the homogeneous information matrix
J_complete. Positive base residual/innovation does not make a homogeneous
error direction observable when its linearized measurement rows vanish.

Consequently the theorem
`q=0 + Delta_g>=theta_E + geomagnetic cutoff + MAGNETIC SERVICE
 => Ker J_complete={0}`
is **not proved** and does not follow from the existing assumptions without a
physical-to-nominal mean constraint.

The neighborhood claim likewise cannot be obtained by continuity. Continuity
would give a uniform gamma_* only after the q=0 matrix is uniformly positive
definite; that premise is exactly what remains missing.

What the new geomagnetic premise does prove is narrower but useful: every
q=0 pure-tilt null word must carry the gravity-sized nominal-AW condition
ZB-1/GF-3. Thus any counterexample is confined to that explicit base-mean
manifold. To exclude it one must use the literal mean adaptation dynamics
over successive words; it cannot be excluded from a single-word homogeneous
information matrix.

## Multi-word mean-dynamic invariance equation for the q=0 manifold

Let `P_B` be the fixed world projection perpendicular to the near-constant
geomagnetic direction and define the compatibility-manifold coordinate
`c_k=P_B(a_w,k-g)` immediately after an accepted accelerometer update.
Persistent q=0 compatibility requires `c_k=0` at every such epoch (up to the
explicit eps_B/lever defect tube).

Propagate from one accelerometer epoch to the next using the literal shipping
order. Write `F_k` for all predictions and any scheduled S pseudo-update
before the next accelerometer correction, and let `L_k` select its world-AW
component. Then before the accelerometer update
`c_(k+1)^-=P_B[L_k F_k x_k-g]`.
The accepted accelerometer update gives
`c_(k+1)^+=c_(k+1)^-+P_B L_k K_a,k r_a,k` (with the literal reset/frame
transport included in L_k). Exact manifold invariance is therefore

`P_B L_k K_a,k r_a,k = -c_(k+1)^-`.                       (MI-1)

Substitute the exact physical residual
`r_a,k=f_phys,k-fhat_k-bhat_a,k`. On the manifold the nominal force has
`P_B(a_w-g)=0`, so MI-1 becomes a linear constraint on the actual varying
physical acceleration plus the endogenous estimator states:

`P_B L_k K_a,k R_wb,k P_B a_phys,k^W
 = -c_(k+1)^- + known BA/parallel/lever terms`.             (MI-2)

Thus an infinite q=0 compatible execution is equivalent to satisfying MI-2 at
every accepted accelerometer epoch together with the S/BA/attitude mean
recursions.

The important point is dimensional. MI-2 is a two-component transverse
constraint on the three-component physical acceleration sample. MARINE MOTION
currently imposes bounds, jerk regularity, bounded velocity/displacement/
potential, and gravity-direction attitude span, but it does not prescribe the
two transverse translational-acceleration components. Therefore MI-2 can in
principle select those components while leaving one physical acceleration
degree of freedom.

The time-varying attitude requirement does not change this count: it changes
the coefficients in MI-2 and the body measurement, but physical translation
remains an independent admissible history. Without a lower excitation
condition on translational acceleration relative to this two-dimensional
constraint, the existing MARINE MOTION class does not algebraically forbid
solutions of MI-2.

Hence the desired theorem
`no infinite execution on P_B(a_w-g)=0`
cannot be derived from the current assumptions by the mean recursion alone.
The recursion instead gives the exact construction equation MI-2 for a
candidate persistent execution. To exclude it one needs either (a) prove that
solutions of MI-2 necessarily violate the all-time velocity/displacement/
potential constraints, or (b) add a physical excitation condition that makes
the transverse physical acceleration incompatible with MI-2 on recurring
windows.

The earlier long-horizon analysis shows MI-2 forcing may be oscillatory and
zero-mean, so (a) is not presently established. Thus the remaining gap is now
an explicit physical-admissibility question for the MI-2-selected acceleration
history, not an estimator observability question.

## Physical admissibility of the MI-2-selected acceleration law

Write MI-2 as
`M_k a_phys,perp,k = h_k + N_k a_phys,parallel,k`,          (PA-1)
where M_k is the literal 2x2 transverse authority obtained from
`P_B L_k K_a,k R_wb,k P_B`, h_k contains the predicted departure from the
compatibility manifold plus BA/lever terms, and the remaining longitudinal
physical acceleration is the free scalar input.

Whenever `sigma_min(M_k)>=d_0>0`,
`a_phys,perp,k=M_k^-1 h_k+M_k^-1 N_k a_phys,parallel,k`.   (PA-2)
Hence the required transverse acceleration is finite and obeys the amplitude
bound
`|a_phys,perp,k| <= d_0^-1(|h_k|+|N_k||a_parallel,k|)`.
On a compact retained history all coefficients are bounded, so MI-2 does not
force acceleration blow-up while transverse authority stays nonsingular.

For jerk, subtract PA-2 at consecutive epochs. With bounded derivatives/
increments of the gains, attitude, BA, lever terms and compatibility departure,
`|Delta a_perp| <= C_M |Delta M| + C_h |Delta h|
                    +C_N |Delta N|+C_u |Delta a_parallel|`. (PA-3)
Thus the physical jerk limit can be met on a sufficiently regular constrained
trajectory; there is no structural lower bound forcing jerk above J_max.
Conversely current assumptions do not guarantee such regularity globally.

For mean and integrated motion, PA-2 is an affine time-varying function of the
free scalar a_parallel. Over a window,
`bar a_perp = bar(M^-1 h)+bar(M^-1 N a_parallel)`.          (PA-4)
One scalar free input cannot generically prescribe both components of the
transverse mean on a single finite window, so zero mean acceleration is not
automatic. But MARINE MOTION does not require zero acceleration mean on every
window; it requires all-time bounded velocity/displacement/potential. Those
conditions require cancellation over longer histories, and a time-varying
scalar input can in principle shape the accumulated vector through the rotating
coefficient `M^-1 N`.

Therefore MI-2 does not structurally violate the current amplitude, jerk, or
finite-window kinematic conditions. The only possible all-time obstruction is
whether its forced transverse component has a nonzero secular mean that cannot
be canceled by the rotating longitudinal channel while respecting the bounds.

Define the long-horizon forced mean
`A_F(T)=integral_0^T M(t)^-1 h(t) dt`
and controllable mean direction
`A_U(T)=integral_0^T M(t)^-1 N(t) u_parallel(t) dt`.
Bounded velocity requires `A_F(T)+A_U(T)=O(1)` as T grows; bounded
displacement/potential impose the corresponding first/second moment
conditions. The current MARINE MOTION and shipping invariants provide no
source-uniform separation of A_F from the reachable set of A_U.

Hence the current assumptions do not exclude the MI-2-selected physical
trajectory. Nor does this prove one exists globally: existence requires a
bounded longitudinal control satisfying the hierarchy of moment conditions
while D_perp remains nonsingular.

This identifies the exact remaining physical theorem: characterize the
long-horizon reachable moment set of the scalar longitudinal channel
`M^-1N`. If the forced moments lie outside it uniformly, O2 closes; if they
lie inside for one forward-complete strict-service history, that history is
the admissible persistent-zero counterexample.

## Fixed-world gravity/geomagnetic mean mismatch

Use the LOCAL GRAVITY reference `g_0^W` and MAGNETIC SERVICE reference
`b_0^W`, and set
`d_0=P_(b_0)^perp g_0^W`.
Then
`||d_0||>=g_min cos(80 deg)`.

For an infinite q->0 compatible execution, the exact accelerometer geometry
implies at every applied epoch
`P_(b_0)^perp a_w,k^W = d_0 + delta_k`,                   (FG-1)
where
`||delta_k|| <= C_g eps_g+C_B eps_B+C_L`.
If the total defect is smaller than `||d_0||`, the nominal AW therefore has
a persistent nonzero fixed-world mean component.

By contrast, physical acceleration satisfies
`integral_0^T a_phys^W dt = v^W(T)-v^W(0)`.
The MARINE MOTION bounded-velocity continuation gives
`(1/T) integral_0^T a_phys^W dt -> 0`.                    (FG-2)
Bounded displacement supplies the corresponding zero long-time mean velocity,
and the bounded-potential condition prevents secular displacement hidden by
finite-window reanchoring.

Hence every persistent compatibility counterexample must realize the strict
mean mismatch
`mean P_(b_0)^perp a_w^W = d_0+O(eps_g+eps_B+C_L) != 0`,
`mean a_phys^W = 0`.                                      (FG-3)

Substituting the exact residual identity into the AW mean recursion shows that
this mismatch must be sustained entirely by the long-time mean of Kalman
measurement corrections:
`mean[Delta a_w,acc + Delta a_w,S + Delta a_w,mag/reset]
 = mean OU loss of d_0 + defects`.                          (FG-4)
The physical acceleration itself cannot supply a DC term.

FG-4 is now the controlling mean-dynamic identity. To exclude infinite
compatibility it remains to prove that, on the q->0 manifold, the literal
measurement-correction sum has zero (or insufficient) fixed-world DC after
all bounded estimator states are telescoped. If that identity holds, the
positive OU loss of d_0 gives an immediate contradiction. If not, the
nonzero correction term identifies the exact mechanism supporting the
persistent invariant zero.

## Long-time telescoping of the literal AW correction supply

Let `Delta_a^meas(k)` denote the sum of all mean increments to the world AW
state at step k caused by actual measurement/pseudo-measurement corrections
(S, accelerometer, magnetic through attitude/reset coupling where applicable),
with prediction treated separately. The exact AW recursion is
`a_(k+1)=phi_k a_k+Delta_a^meas(k)`.                       (TEL-1)

Summing,
`a_N-a_0=-sum_(k<N)(1-phi_k)a_k+sum_(k<N)Delta_a^meas(k)`. (TEL-2)
For bounded a_w, division by elapsed time makes the left side vanish. On the
persistent q->0 manifold FG-1, projection on d0 gives a strictly positive OU
loss rate (minus the explicit field/gravity/lever defects). Therefore
`mean_T P_B Delta_a^meas = mean_T (1-phi) P_B a_w >0`.     (TEL-3)

So the total measurement-induced AW correction does **not** telescope to zero;
bounded AW proves the opposite: it must have a positive mean exactly balancing
OU decay.

Can bounded v,p,S,BA,attitude force the individual correction sum to zero?
No. Their exact telescoping equations constrain different gain-weighted
combinations of the same innovations. For example the S pseudo-update gives
increments `K_awS(-S)`, `K_vS(-S)`, `K_pS(-S)`, etc.; bounded S constrains
the S-state row of the total update, not the AW row. Cross-covariance gains
provide no identity equating their long-time means. Likewise bounded BA and
attitude constrain their own gain rows. Joseph covariance updates do not
supply a mean conservation law.

The accelerometer correction also need not have zero mean when physical
acceleration has zero mean, because
`r_a=f_phys-fhat-bhat_a`; persistent compatibility gives fhat a fixed
gravity-related offset. Thus zero mean physical acceleration is compatible
with a nonzero mean residual and hence a nonzero mean AW correction.

Therefore the hoped-for contradiction
`bounded states => mean measurement AW correction=0`
is false. TEL-3 is a balance law, not an impossibility theorem.

The exact remaining requirement for a persistent zero is that the gain-weighted
correction supply equal the OU loss while all other state-row balance equations
hold simultaneously. Stack the bounded-state telescoping equations:
`0 = -Lambda xbar + Gbar_r + Gbar_S + defects`,            (TEL-4)
where each row uses its literal gain history. Persistence is feasible iff this
coupled mean-balance system has a solution consistent with the zero-mean
physical acceleration and compatibility constraints.

No existing Kalman identity makes TEL-4 inconsistent. Thus LOCAL GRAVITY and
near-constant geomagnetism expose a large required correction supply, but
bounded estimator states alone do not eliminate it. A proof of impossibility
must establish a rank/range separation for the **stacked mean correction
matrix** in TEL-4; otherwise TEL-4 is the algebraic mechanism supporting the
persistent invariant zero.

## Stacked long-time mean correction range test

Stack the bounded mean-state rows that participate in the compatibility
balance, `x_m=(a_w,v,p,S,b_a,theta,...)`. At each correction epoch the same
innovation source enters all rows through one gain matrix. Over a long horizon
define the gain-weighted source averages without separating linked histories:
`z_a=(1/T)sum K_a,k r_a,k`, `z_S=(1/T)sum K_S,k(-S_k)`, and analogous
magnetic/reset contributions. The bounded-state telescoping equations have
the block form
`b_req = R_a zeta_a + R_S zeta_S + R_m zeta_m`,            (RG-1)
where `b_req` contains the deterministic prediction drifts. Its AW component
contains the positive OU replenishment of d0; its LIN components contain the
mean deterministic chain couplings.

The key point is that the source variables are not arbitrary vectors:
`zeta_a` is generated by one physical residual history and `zeta_S` by the
endogenous S history. Nevertheless the literal Kalman gains are full NX-by-3
matrices. Algebraically, the stacked correction range contains the columns of
`[K_a,K_S,K_m]` accumulated over time.

There is no structural row identity in the shipping filter making the AW row
a linear combination of the v,p,S rows. The covariance cross blocks are
independent state covariances subject only to PSD/Joseph recursion. Therefore
a nonzero AW mean correction can coexist algebraically with zero secular
v,p,S balance: S and accelerometer sources provide different gain columns and
can cancel LIN rows while adding in the AW row.

Equivalently, no nonzero left annihilator l is known (or implied by the state
model) such that
`l' [R_a R_S R_m]=0`
for every admissible history while
`l' b_req !=0` due to the geomagnetic d0 term. Without such an annihilator,
the Fredholm/range test does not exclude the required correction supply.

Indeed the four-S complete-word injectivity result points in the opposite
direction for the finite-horizon linear system: the S-chain source/action
columns are sufficiently independent to determine LIN/AW roots rather than
imposing a conserved mean quantity. Measurement corrections break the open
chain's simple integral conservation.

Thus the stacked **algebraic** mean correction matrix does not yield a
source-uniform contradiction from the current shipping structure. The
required AW DC vector is not structurally outside its range.

What remains is the dynamic/source constraint: whether the particular source
averages needed by RG-1 can arise from one same-history residual/S sequence
whose physical acceleration has zero long-time mean and satisfies MARINE
MOTION. That is not a static rank question. Any proof must use temporal
correlations/moments of the source sequence; static range has enough degrees
of freedom.

Therefore the rank/range test is resolved negatively as an exclusion route:
there is no shipping left-null invariant presently separating the required
geomagnetic AW replenishment from the correction range. Persistent-zero
exclusion, if true, must come from temporal moment constraints rather than
instantaneous/mean correction rank.

## Temporal-moment constraints from bounded physical motion

For one world component, physical kinematics satisfy
`a=dv/dt`, `v=dp/dt`. Bounded velocity gives the zeroth acceleration
moment
`A_0(T)=integral_0^T a(t)dt=v(T)-v(0)=O(1)`, so `A_0(T)/T->0`.

Integration by parts gives the first acceleration moment
`A_1(T)=integral_0^T (T-t)a(t)dt
       =p(T)-p(0)-T v(0)`.                                 (TM-1)
Thus bounded displacement does **not** make A_1 bounded unless the initial
velocity is zero; it fixes its secular part exactly. Equivalently, after
subtracting the boundary velocity term, the centered first moment is O(1).

The proof's stronger bounded-potential condition
`|integral_(t1)^(t2) p(t)dt|<=P_AC`
controls the next integrated displacement moment and prevents repeated
reanchoring/secular offsets. These are boundary/moment identities, not a rule
that every finite-window acceleration mean is zero.

Apply the same moment operators to the MI-2-selected acceleration law
`M a_perp=h+N a_parallel`. The zeroth and centered-first moment constraints
become
`integral M^-1 h + integral M^-1 N a_parallel = O(1)`,     (TM-2)
and
`integral (T-t)[M^-1 h+M^-1N a_parallel]
   = prescribed boundary term + O(1)`.                     (TM-3)
The bounded-potential condition adds the corresponding second integrated
constraint.

Because the control direction `M^-1N` is time varying under MARINE MOTION,
one scalar function a_parallel(t) has infinitely many temporal degrees of
freedom and can, in principle, satisfy finitely many vector moment constraints
over successive windows. There is no dimensional contradiction analogous to
a single scalar constant trying to cancel two fixed vectors.

Most importantly, the fixed nominal AW offset d0 does not appear directly as
a physical acceleration DC term in TM-2/TM-3; it enters h through the Kalman
correction requirement. The static range test already showed those correction
columns can algebraically supply that offset. Time weighting alone does not
create a sign-definite invariant.

Therefore bounded displacement and bounded potential strengthen the admissible
source constraints but still do not, from the current identities, force an
uncancelable secular moment. A contradiction would require a source-uniform
moment-separation theorem for the rotating scalar channel `M^-1N` (for
example, a nonzero left functional annihilating all admissible control moments
but not the forced h moments). No such functional follows from the current
MARINE MOTION assumptions.

This resolves the naive first-moment route: it supplies exact boundary
conditions but not yet O2 exclusion. The remaining question is a temporal
controllability/moment problem for the MI-2 law, not a missing integration-by-
parts identity.

## Temporal controllability of the rotating scalar MI-2 channel

Set
`f(t)=M(t)^-1 h(t)` and `g(t)=M(t)^-1 N(t) in R^2`, so the selected
transverse physical acceleration is
`a_perp(t)=f(t)+g(t)u(t)`.                                 (TC-1)

Consider a linear time-local functional
`L_w[a]=integral w(t)'a(t)dt`.
For L_w to annihilate the contribution of **every** admissible scalar control
u with compact support inside a regular interval, the fundamental lemma of
the calculus of variations requires
`w(t)'g(t)=0` almost everywhere.                            (TC-2)
In two transverse dimensions, where g(t) is nonzero, every such annihilator is
`w(t)=lambda(t) J g(t)`, with J the 90-degree rotation.

Therefore the uncontrollable instantaneous component of the MI-2 law is
exactly
`f_perp_ctrl(t)=P_(g(t))^perp f(t)`, or equivalently the scalar
`chi(t)=(Jg(t))' f(t)`.                                   (TC-3)

This resolves the abstract static/moment controllability question: arbitrary
time variation of the scalar longitudinal input can generate only the
pointwise line span{g(t)}; it can never alter chi(t). Time weighting and
higher moments do not change this local annihilator.

A source-uniform exclusion theorem would follow if persistent compatibility
forced a sign-definite/nonzero long-horizon moment of chi, for example
`liminf_(T->infinity) |integral_0^T chi(t)dt|/T >= chi_0>0`,
because bounded physical velocity requires every fixed-world acceleration
component to have zero long-time mean.

However chi is expressed through the same-history Kalman gains/cross
covariances in M,N,h. The fixed gravity/geomagnetic vector d0 enters h, but
nothing presently proved fixes the sign of its projection on Jg(t).
MARINE MOTION gravity span rotates the coefficients and can make chi change
sign. Hence no positive chi_0 follows from the current assumptions.

Conversely, if `chi(t)=0` identically along a compatible history, then the
forced transverse acceleration lies pointwise in the controllable line and a
scalar u(t) can cancel it exactly at the acceleration level. More generally,
zero long-time moments of chi remove the velocity obstruction, though
displacement/potential still impose integrated conditions.

Thus the temporal controllability problem has a precise invariant:
`chi(t)=det[g(t),f(t)]`.                                   (TC-4)
The current theorem assumptions neither bound chi away from zero nor force a
nonzero secular moment. This is the remaining physical-to-estimator linkage
needed for O2. Any additional assumption, if ultimately required, should be
stated physically so that it implies a nonzero recurring chi-action; it should
not directly constrain Kalman gains.

## Determinant reduction of the temporal uncontrollable component

For the 2x2 invertible transverse authority M,
`chi=det(M^-1 N,M^-1 h)=det(N,h)/det(M)`.                  (CH-1)
Therefore M cannot create or remove the zero of chi; it only scales and
orients it. The decisive geometry is the pair (N,h) before transverse
inversion.

N is the effect of the one remaining longitudinal physical-acceleration
component on the two compatibility-maintenance equations. h is the forced
departure produced by OU decay, scheduled S pseudo-updates, BA/lever terms,
and the fixed gravity/geomagnetic offset d0.

There is no shipping identity forcing N and h to be nonparallel. In
particular, h contains gain-weighted estimator-state terms, not just d0.
Those terms evolve with the same covariance/attitude history that determines
N. The fixed world vectors g0,b0 constrain one component of h but do not fix
its direction in the two-dimensional compatibility-equation space.

Recurring gravity-direction span changes N and h continuously but likewise
does not prohibit repeated or persistent collinearity. A rotating pair of
vectors may remain parallel for all time. Thus
`Delta_g>=theta_E` plus fixed noncollinear g0,b0 does not imply
`|det(N,h)|>=chi0` or a nonzero signed average.

Consequently the literal shipping structure and current physical assumptions
do **not** force a nonzero recurring chi-action. The possibility
`det(N,h)=0` is an exact codimension-one compatibility condition that the
remaining physical/estimator history can in principle satisfy; no existing
invariant excludes it.

This answers the temporal-controllability fork negatively as an exclusion
route: chi need not be separated from zero by the present assumptions.
A zero-moment chi is likewise not excluded because no sign law is available.

The pathological mechanism is now explicit: a correlated marine translation
history can choose its longitudinal component so that the forced departure h
lies in the instantaneous longitudinal-control image N, while the two
transverse acceleration components are fixed by MI-2. If this condition is
maintained with bounded moments, the q=0 compatibility manifold can persist.

This is still not a constructed global execution, but it shows that no theorem
based only on gravity span, fixed local gravity, near-constant nonvertical
geomagnetism, and the existing bounded-motion conditions can derive a
source-uniform chi floor without an additional physical restriction linking
translation to attitude. Any such added restriction should target this
correlated-translation degeneracy, not generic acceleration magnitude.

## Literal frozen S-to-S transverse DC calculation

For one fixed transverse world direction and one frozen/slowly varying
coefficient pattern, use the shipping per-axis LIN order
`x=(v,p,S,a)'`. Let `Phi` be the exact analytic 4x4
`PhiAxis4x1_analytic(tau,Delta)` propagation from one scheduled S epoch to
the next, including all intervening prediction time. Let `e_S=(0,0,1,0)'`.

At the scheduled S=0 pseudo-measurement the scalar projected gain column is
`k_S=(k_vS,k_pS,k_SS,k_aS)'`, so
`x^S=(I-k_S e_S') Phi x`.                                 (SDC-1)

Collect the net projected accelerometer mean correction over the same lifted
interval as one scalar innovation/control u with frozen lifted gain column
`k_A=(k_vA,k_pA,k_SA,k_aA)'`. Then the exact frozen lifted mean map is
`x_+ = A_S x + B_S u`,                                    (SDC-2)
`A_S=(I-k_S e_S')Phi`, `B_S=k_A`,
with additional accelerometer epochs represented by the chronological product
and summed input columns; the one-column form is the single-effective-input
case.

A bounded equilibrium satisfies
`(I-A_S)x_*=B_S u_*`.                                     (SDC-3)
Compatibility requires
`e_a' x_*=g_perp`, `e_a=(0,0,0,1)'`, with
`g_perp>=g_min cos80` after the declared defects.

If `I-A_S` is invertible, define the literal frozen DC gain
`G_DC=e_a'(I-A_S)^-1 B_S`.                                (SDC-4)
Then the required accelerometer innovation is exactly
`u_*=g_perp/G_DC`                                         (SDC-5)
provided `G_DC!=0`. The remaining equilibrium states are
`x_*=(I-A_S)^-1B_S g_perp/G_DC`.

Therefore a bounded frozen equilibrium is **not** excluded merely by the OU
integrator/S pseudo-update structure. It exists algebraically whenever
`I-A_S` is invertible and `G_DC!=0`; its physical admissibility is decided
by the magnitude/sign of u_* and by whether the resulting physical
acceleration history satisfies MARINE MOTION. If `G_DC=0`, no finite scalar
accelerometer innovation can sustain nonzero g_perp in this frozen
single-input model. If `I-A_S` is singular, use the left-null/Fredholm
conditions together with the compatibility row.

The Kalman code supplies no structural identity forcing `G_DC=0`.
The S pseudo-update gain k_aS can have either sign through P_aw,S, and the
accelerometer AW gain k_aA is generally nonzero. Thus outcome A (no bounded
equilibrium) cannot be proved from architecture alone; outcome B is
algebraically generic for a frozen regular coefficient set.

For the real MOVING theorem this is only a local diagnostic: gains, attitude
and physical acceleration vary, and MARINE MOTION forbids replacing the
whole execution by constant forcing. But SDC-4--SDC-5 answer the requested DC
question: the literal S adaptation does not inherently reject a gravity-sized
nominal AW equilibrium. The required innovation is the reciprocal lifted DC
gain times g_perp. A persistent counterexample would have to realize the
time-varying analogue of SDC-5 with zero long-time physical acceleration mean.

## Physical realizability of the time-varying required innovation

The lifted control u_j in SDC-5 is an accelerometer **innovation**, not physical
acceleration. On the compatibility manifold its projected value has the exact
form
`u_j = H_j a_phys,j^W - H_j a_w,j^W - H_b bhat_a,j + d_g,j`, (PRDC-1)
with the literal frame/lever terms included in H_j,d_g,j. Because persistent
compatibility pins the relevant nominal AW component near the fixed d0, a
nonzero mean innovation is naturally produced even when physical acceleration
has zero long-time mean.

Let `u_req,j` be the innovation required by the time-varying lifted
compatibility recurrence (the analogue of g_perp/G_DC). Solving PRDC-1 gives
the required physical acceleration component
`H_j a_phys,j^W = u_req,j + H_j a_w,j^W + H_b bhat_a,j-d_g,j`. (PRDC-2)

Therefore zero-mean physical acceleration is feasible iff the long-time mean
of the right side of PRDC-2 vanishes (with the full vector equations and
moment constraints), not iff mean(u_req)=0. The fixed nominal AW term can
cancel a nonzero innovation mean.

Amplitude and jerk feasibility follow locally from bounded/invertible lifted
coefficients as in PA-2/PA-3. Recurring gravity-direction span can be imposed
independently through a smooth periodic/aperiodic attitude history; it changes
the coefficients in PRDC-2 but does not by itself force a nonzero physical
acceleration mean.

Thus outcome B is **physically possible in principle** under the present
assumptions: nothing requires the time-varying innovation needed to sustain
the biased nominal AW to correspond to a DC physical acceleration. A
zero-mean, bounded, jerk-limited physical acceleration can generate a
nonzero-mean innovation because the nominal prediction itself is biased by
the gravity-sized compatibility offset.

This is not yet an explicit global shipping trajectory, because the same
history must also make the covariance/gain sequence generate u_req and satisfy
all higher moment/service conditions. But the proposed physical-realizability
obstruction (zero mean acceleration versus nonzero required innovation) is
not valid.

Consequently current MARINE MOTION does not exclude outcome B at the level of
physical kinematics. A rigorous counterexample now requires solving the
self-consistent covariance/gain recurrence; a proof of O2 would require an
additional invariant showing PRDC-2 cannot be self-consistent, not merely
bounded-motion kinematics.

## Exact chronology test for a D_perp covariance floor

Shipping prediction first propagates the LIN/AW covariance and injects fresh
AW process covariance. Pending AW covariance inflation is then applied.
If due, the S=0 pseudo-measurement performs a Joseph covariance update **before**
the next accelerometer gain is formed. Thus the accelerometer numerator uses
the post-S covariance, not the raw post-prediction covariance.

Let `N_a=P C_a'` denote the accelerometer gain numerator and project its AW
rows onto the transverse compatibility channel. Before S correction the fresh
prediction contributes a positive AW marginal term to N_a. The S covariance
update is
`P^S=P^- - P^-H_S'(H_SP^-H_S'+R_S)^-1 H_SP^-`
(in exact covariance form). Hence
`N_a^S=N_a^- - P^-H_S' S_S^-1 H_S P^- C_a'`.              (DF-1)

The subtractive term has no sign/alignment restriction relative to the
projected fresh AW contribution. PSD only guarantees the **full** covariance
remains PSD. Therefore DF-1 can algebraically cancel a projected AW gain
numerator while all innovation covariances remain SPD.

Subsequent accelerometer gain multiplication by `S_a^-1` cannot restore a
lost projected rank. Thus there is no source-uniform positive D_perp floor
derivable solely from fresh AW process covariance.

Conversely DF-1 does not force cancellation: exact rank loss is an algebraic
surface in the post-prediction covariance/geometry variables. The S update can
move toward or away from it. Prediction at the next step again injects fresh
AW variance. No monotone recurrence drives the system to this surface.

Therefore the actual evolving covariance recursion neither forces
`D_perp->0` nor supplies a positive invariant floor by operation-wise
structure. The rank-loss set is a reachable-looking algebraic boundary, not an
attractor or repeller established by the current invariants.

This means the requested global yes/no cannot be decided from covariance
recursion identities alone. To prove a persistent counterexample one must show
one compatible trajectory avoids the DF-1 cancellation surface for all time;
to prove O2 one must show every compatible trajectory hits it (or another
admissibility boundary). Neither statement follows from shipping covariance
algebra.

## Signed determinant barrier test for the S-induced cancellation surface

Let `d(P,x)=det D_perp(P,x)` on a regular fixed-event compatibility branch.
The cancellation surface is `Z={d=0}`. A compact one-sided invariant region
would require a source-uniform sign condition for the constrained update of d
at Z.

The post-S covariance is the smooth map
`P^S=P-PH_S'(H_SPH_S'+R_S)^-1H_SP`.
On the SPD innovation domain its differential with respect to P is finite and
generically nonzero. D_perp also depends smoothly on the nominal geometry and
the subsequent transport. Therefore d is a smooth scalar function on each
regular branch.

At a regular rank-one point of Z, the first variation is
`delta d = tr(adj(D_perp) delta D_perp)`.                  (BAR-1)
The admissible compatible history supplies variations through preceding
attitude/geometry, prediction covariance and the remaining physical input.
No shipping identity constrains BAR-1 to one sign. In particular the S
subtraction term and fresh AW prediction term enter delta D_perp with opposite
possible projected orientations. Hence Z is not proved invariant and no
one-sided inward barrier follows from the covariance algebra.

But transversality of Z is not enough to construct an infinite avoiding
trajectory. A codimension-one surface can be crossed by generic trajectories;
local freedom can choose a side only while the compatibility control remains
regular. To prove eternal avoidance one needs a controlled-invariance theorem
for a compact subset `|d|>=d0`, including physical moment/service bounds.
No such inward condition is supplied by BAR-1.

Thus the proposed signed-determinant barrier route does not close either
branch. It establishes that exact S-induced cancellation is generically a
crossable hypersurface rather than an estimator-internal invariant, which
makes mandatory rank loss implausible, but it does not prove a forward-
complete nonsingular compatible trajectory.

The mathematically honest status is therefore: current assumptions do not
force every compatible trajectory to hit Z, and local compatible trajectories
can avoid/cross Z; global eternal avoidance remains an existence problem.
A proof of it requires an additional recurrence/compactness mechanism beyond
the local determinant sign.

## Periodic compatible physical/mean construction: closure conditions

Choose a smooth periodic physical attitude R(t) with period T_p and gravity
span strictly above theta_E on every complete T_E window (for example a
nondegenerate rocking motion with period/cadence chosen so the rolling-window
span condition holds). Keep local gravity and geomagnetic world references
fixed within their declared envelopes. This gives strict recurring magnetic
geometry away from the dip-pole exclusion.

Choose physical translation periodic with the same T_p. Periodicity of
position and velocity is equivalent to the acceleration waveform satisfying
the zeroth and first moment closure conditions over one period:
`integral_0^Tp a(t)dt=0`,
`integral_0^Tp (Tp-t)a(t)dt=0`.                            (PC-1)
Then velocity/displacement repeat and the bounded-potential condition holds
for the periodic zero-mean displacement after choosing its reference.

Now impose q=0 compatibility at every accelerometer epoch. For a given
covariance/gain history this determines the two transverse acceleration
components through MI-2, leaving one longitudinal scalar sequence u_j. Over
one physical period, PC-1 gives vector closure constraints on u_j. Because
there are many accelerometer epochs per period, u_j supplies many scalar
degrees of freedom; the closure problem is a finite linear/nonlinear
boundary-value system rather than an overdetermined single-control equation.

However the required transverse waveform depends on P through the gains.
Therefore a periodic physical/mean orbit cannot be constructed independently
of covariance and then repeated while P changes: changing P changes MI-2 and
hence the acceleration needed for exact compatibility. This invalidates the
shortcut of fixing a periodic physical/mean period while allowing arbitrary
nonperiodic covariance evolution.

A genuine repeated exact-compatible period requires either:
(a) the relevant gain/authority sequence repeats (not necessarily full P), or
(b) the physical acceleration waveform is adapted from period to period to the
evolving P, in which case the physical/mean orbit is not periodic.

Thus continuity of the finite-period covariance map alone cannot preserve a
fixed determinant margin under repetition. The proposed construction reduces
back to a joint recurrence of the **relevant gain projection** and the
periodic closure constraints. Full covariance recurrence is stronger than
necessary, but recurrence of the gain quantities entering MI-2 is necessary
for an exactly repeated physical/mean orbit.

This identifies the minimal periodic fixed-point variables: nominal compatible
mean/tuner state plus the projected covariance blocks entering K_a,K_S and
D_perp, not the entire P. A periodic counterexample may be sought on that
reduced map. No theorem currently proves that reduced map has a fixed point.

## Parametric compatible-period closure and Riccati dichotomy

Fix a smooth periodic attitude/service/event template with strict MOVING and
magnetic-service margins. Over one period let P0 be the entering covariance,
z the entering compatible nominal mean/tuner state, and u the finite vector of
free longitudinal physical-acceleration controls at the accepted
accelerometer epochs (with smooth interpolation between samples).

Collect all one-period constraints in
`F(P0,z,u)=0`: exact compatibility at every accelerometer epoch, terminal
mean/tuner closure, zero velocity/displacement period moments for physical
translation, and any scheduler phase closure. The two transverse physical
acceleration components are eliminated by MI-2 wherever D_perp is nonsingular.

If at one regular solution `(Pbar,zbar,ubar)` the Jacobian
`D_(z,u)F` is onto, the implicit-function theorem gives a local smooth
selection `z(P),u(P)` for every P near Pbar. Strict physical/service margins
and D_perp nonsingularity persist by continuity. Thus compatible mean/physical
period closure can be maintained while P varies locally; the waveform is
allowed to adjust from period to period.

Define the exact period covariance map along that selected compatible period
by
`P_+=R(P)`.                                                (PR-1)
Iterate `P_(n+1)=R(P_n)`, re-solving z(P_n),u(P_n) each period. This yields
a forward compatible execution as long as P_n remains in the IFT domain (or
overlapping continuation charts) and no physical/service/rank boundary is
reached.

There is then a rigorous dichotomy **conditional on global continuation of
the closure charts**:
(A) if {P_n} is unbounded while the compatible physical/mean execution remains
admissible, the desired uniform A21 covariance ceiling is false directly;
(B) if {P_n} is bounded, compactness gives accumulation points. A fixed or
recurrent covariance orbit is not automatic from boundedness alone, but the
omega-limit set is a nonempty compact invariant set for a continuous R on a
closed continuation domain. D_perp remains separated from zero if the
continuation domain was constructed with that margin.

This formulation avoids requiring an identical physical waveform or a
periodic P. It also exposes the one missing existence premise: a **single
regular seed solution** with onto D_(z,u)F and D_perp!=0. Constructor
covariance is not enough unless it lies on an actual recurring A21-compatible
period seed.

Therefore the bounded/unbounded Riccati dichotomy is ready once a regular
parametric period seed is proved. The remaining construction task is finite:
exhibit one admissible period and show the boundary-value Jacobian with respect
to the free mean/control variables has full row rank.

## Regular compatible-period seed: Jacobian structure and remaining existence gap

Choose a smooth periodic rocking attitude template with strict
`Delta_g>theta_E`, strict MAGNETIC SERVICE, and scheduler phases repeating
after T_p. Parameterize physical acceleration at the accepted accelerometer
epochs by three components per epoch and use a smooth periodic interpolation.

Before eliminating controls, the boundary-value equations consist of:
(i) two compatibility equations per accelerometer epoch;
(ii) physical velocity/displacement closure over the period;
(iii) nominal mean/tuner/scheduler closure.

At a point with D_perp nonsingular, the derivative of block (i) with respect
to the two transverse acceleration components at each epoch is block lower
triangular in chronology with invertible diagonal blocks D_perp,k. Hence all
epoch-wise compatibility equations are locally solvable and can be eliminated.

After this elimination, many longitudinal acceleration samples remain. Their
derivatives with respect to terminal velocity and displacement are the usual
discrete moment rows `sum w_k u_k` and `sum (T_p-t_k)w_k u_k`. With at
least two distinct control epochs and nonzero longitudinal direction these
rows are independent for the scalar longitudinal component; rotating attitude
can provide vector closure authority over a full period. Thus kinematic
closure is generically controllable, but a source-uniform full-rank proof
requires an explicit lower bound on the corresponding sampled moment matrix.

The nominal mean closure derivative with respect to the initial mean is
`Phi_mean(T_p)-I` after the compatibility controls are eliminated. This block
is not generically invertible: neutral attitude/reference coordinates and
integrator coordinates can have unit Floquet multipliers. Closure may use
remaining control degrees of freedom, but surjectivity then depends on the
full endpoint controllability matrix of the reduced mean recursion.

Tuner closure is more problematic: tuner states are deterministic functions
of the measurement/mean history and may include clamps/discrete branches.
They are not free initial coordinates on a fixed branch unless the branch is
inside a smooth unclamped region. No current theorem proves a periodic tuner
orbit for the proposed compatible forcing.

Therefore the triangular Jacobian argument proves only the first block
(compatibility) rigorously from D_perp!=0. It does not yet prove full
D_(z,u)F surjectivity or exhibit a seed. Claiming a regular seed from dimension
count would be incorrect.

The exact finite obligation is now:
(a) choose a concrete smooth tuner/event branch;
(b) prove full-rank endpoint controllability of the reduced longitudinal-input
mean recursion over one period, including v,p,S,AW/BA/attitude closure;
(c) show the tuner map has a fixed point on that branch;
(d) select an entering covariance with D_perp!=0.

Until (a)-(d) are supplied, the regular compatible-period seed remains
unproved.

## Endpoint controllability: correct reduced target

After imposing the periodic physical attitude/service template and eliminating
the two transverse physical acceleration components by compatibility, do not
ask the remaining scalar input to control every estimator coordinate. Attitude
is prescribed/closed by the compatibility+magnetic construction, scheduler
phase is chosen periodic, and fixed filter parameters need no control.

Let y collect only the residual mean closure coordinates not already fixed by
those constraints (LIN/S and any active BA scalar combinations). The reduced
linearized endpoint map is
`delta y_N = Phi_y delta y_0 + C_N delta u`,
`C_N=[Phi(N,k+1)b_k]_(k in acc epochs)`,                  (ECM-1)
where b_k is the effective column after transverse compatibility elimination.

Full endpoint controllability is exactly
`rank C_N = dim(y_unfixed)`.                              (ECM-2)

The open OU/LIN chain driven by an acceleration input has a confluent
Vandermonde/Hermite controllability structure: distinct input epochs generate
independent endpoint moments in a_w,v,p,S. Four distinct effective epochs are
enough per scalar chain when the direct AW component of b_k is nonzero. This
is the same {1,t,t^2,psi_tau(t)} structure used by the four-S injectivity
lemma, now transposed as a controllability statement. Thus the LIN/AW
subsystem is controllable on a regular branch with nonzero effective AW input.

Active BA/attitude closure is not supplied by that scalar chain. Those
coordinates are already constrained by exact compatibility and the prescribed
periodic attitude; requiring them again in y double-counts constraints. The
remaining BA amplitude on the compatibility line decays and is zero in the
q=0 periodic construction. Therefore the legitimate q=0 seed closure target
is the LIN/AW chain, for which ECM-2 holds with four distinct regular
accelerometer epochs.

Hence endpoint controllability of the **correct reduced q=0 mean closure**
closes conditionally on the effective longitudinal input having a nonzero AW
component at four distinct epochs. That nonzero component is the same regular
gain condition underlying D_perp and is open.

## Unified finite-error physical-to-nominal vector-diversity lemma

Define on a source-qualified complete history window W the literal chronological
joint action J_vec(W,e,n): applied accelerometer and magnetic innovation losses
plus the actual AW/BA process action and S pseudo-measurement action, with
resets/transports and same-history covariance whitening retained. Here e is
the finite attitude error and n denotes admissible estimator nuisance
trajectories (AW, BA, LIN/S and associated mean corrections).

The correct theorem is a separation statement, not unconditional coercivity:

For every retained radius r and every epsilon>0, there exists
`c_vec(epsilon,r,theta_E)>0` such that
`J_vec >= c_vec`
whenever the finite-error state is at physical-metric distance at least epsilon
from the exact compatibility set C.                            (VD-1)

Equivalently, any normalized sequence with J_vec->0 has, after same-history
compactness extraction, a limiting trajectory in C.           (VD-2)

The proof is by contradiction. Physical histories are compact from the
MARINE MOTION acceleration/jerk, velocity, displacement/potential and
gravity-span assumptions; local gravity and geomagnetic histories converge to
their declared compact field classes; finite estimator/tuner/covariance
variables are retained. Zero limiting magnetic loss enforces magnetic-axis
attitude compatibility. Zero accelerometer loss, together with zero AW/BA/S
process/pseudo action, gives the previously classified tilt/BA compatibility
relation. Four distinct S observations remove independent LIN/AW homogeneous
roots. Hence the only zero-action limit is C. Distance>=epsilon is closed and
disjoint from C, contradiction.

This lemma is valid only if event-boundary lower semicontinuity is carried with
the raw/variational action; do not use pseudoinverse constant-rank strata.

**Capture use.** If every attitude error with tilt>=7 degrees is separated by
a positive physical-metric distance from C on the H18/capture branch, VD-1
gives positive recurring action and finite entry into the refinement tube.
That separation is an additional branch-specific fact and must be proved; C
must not contain a >=7-degree captured-domain alias.

**O2 use.** As epsilon->0, VD-1 does not yield a positive floor on C. Instead
VD-2 identifies every vanishing-action sequence with the explicit
compatibility manifold. The remaining O2 question is then exactly whether the
closed-loop mean dynamics can remain in C indefinitely. Thus the lemma
unifies the reduction but does not by itself prove escape from C.

This is the strongest common theorem supported by the existing nullspace and
compactness machinery. Claiming that its small-error limit automatically
closes O2 would be false; the compatibility manifold is deliberately retained.

## H18 finite-error compatibility outside the 7-degree tube

On H18 the accelerometer-bias estimate is held (and its decoupled covariance
block is unchanged on feasible held segments). Zero magnetic loss between true
and nominal records with fixed/near-fixed world b0 restricts the attitude
discrepancy to a rotation `Q_delta` about b0.

For zero accelerometer loss at every applied epoch, the same held body-frame
BA correction must satisfy the difference between the true and nominal
specific-force records. Ignoring only the declared bounded field/lever defects,
the required held BA is
`b_H = R_true(t)[a_phys(t)-g0]
       -R_hat(t)[a_w_hat(t)-g0]`.                           (H7-1)
With `R_hat=R_true Q_delta` in a consistent world/body convention, this
becomes a time-varying function of R_true(t), physical acceleration and
nominal AW unless Q_delta=I or the histories satisfy a special compatibility
relation.

A constant held BA therefore cannot generically absorb a nonzero magnetic-axis
rotation over a moving history. However the H18 nominal AW state is not held:
it is an estimator nuisance trajectory and can vary through prediction and
measurement corrections. Solving H7-1 for nominal AW gives
`a_w_hat(t)-g0 = Q_delta^-1 [a_phys(t)-g0-R_true(t)^-1 b_H]`. (H7-2)
Thus for any fixed Q_delta and b_H there is algebraically a time-varying
nominal AW trajectory that makes the accelerometer residual zero.

The H18 compatibility question therefore reduces to whether the actual H18
OU/LIN/S mean dynamics can realize H7-2. The held-BA property alone does not
exclude rotations >=7 degrees. Fixed g0,b0 and gravity-direction span make the
required AW trajectory nontrivial, but they do not bound its amplitude/action
away from the allowed H18 nuisance class without using the AW/S dynamics.

Consequently
`C_H18 intersect {tilt>=7deg}=empty`
is **not proved** by magnetic geometry plus held BA alone. The same
physical-to-nominal AW bridge remains. A finite-error capture proof must include
the H18 AW/S process action in the variational lemma; if zero joint action is
assumed, then H7-2 must also satisfy zero AW process/S action, which is much
more restrictive and is the correct next test.

## H18 zero-action OU substitution

For an exact zero-loss H18 alias let Q be the constant attitude discrepancy
rotation about the fixed world magnetic direction b0, and let b_H be the held
accelerometer-bias estimate in the body/error convention used by H7-2. The
accelerometer zero-loss identity is
`ahat_w-g0 = Q^-1[a_phys-g0-R_true' b_H]`.                 (HO-1)

Zero AW process action requires the homogeneous OU mean law. On an interval
with constant tau,
`d ahat_w/dt = -ahat_w/tau`.                               (HO-2)
Substitution and differentiation give the necessary physical law
`dot a_phys - d(R_true' b_H)/dt
 = -(1/tau)[a_phys-g0-R_true'b_H+Q g0]`.                    (HO-3)
For piecewise-varying shipping tau the same relation holds intervalwise with
the corresponding tau(t).

More decisive is the long-time mean of HO-1+HO-2. Any bounded homogeneous OU
trajectory has time average zero:
`mean ahat_w =0`.
Bounded physical velocity gives
`mean a_phys=0`.
Hence an infinite zero-action alias would require
`0-g0 = Q^-1[0-g0-mean(R_true' b_H)]`, or
`mean(R_true' b_H)=(Q-I)g0`.                              (HO-4)

The left side has norm at most |b_H|. Therefore a necessary condition is
`|(Q-I)g0| <= |b_H|`.
For a rotation angle alpha about b0,
`|(Q-I)g0|=2 |P_bperp g0| sin(|alpha|/2)`.
Using |P_bperp g0|>=g_min cos80deg gives
`2 g_min cos80deg sin(|alpha|/2) <= |b_H|`.               (HO-5)

Thus every H18 zero-action magnetic-axis alias obeys the explicit angle bound
`|alpha| <= 2 asin(|b_H|/[2 g_min cos80deg])`, provided the
argument is <=1.                                                    (HO-6)

Using the declared physical/held accelerometer-bias bound B_H gives the
source-uniform ceiling
`alpha_H,max = 2 asin(B_H/[2 g_min cos80deg])`.            (HO-7)

With the proof's recorded physical accelerometer-bias bound
B_H=0.22516660498395405 m/s^2 and g near standard gravity, the denominator
2*g*cos80deg is about 3.406 m/s^2, giving alpha_H,max about 0.132 rad,
approximately 7.58 degrees. Therefore the current coarse bias/nondip-pole
bounds do **not** yet exclude every 7-degree alias; they miss by roughly
0.58 degree before eps_g/eps_B/lever defects.

This is nevertheless a sharp quantitative result. To exclude >=7 degrees one
needs either the actual tighter held-BA bound on H18, a stronger inclination
domain than 80 degrees, or use the attitude-span average in HO-4 to improve
`|mean(R' b_H)|<=|b_H|` strictly. The last option uses existing MOVING
excitation and is preferable to changing assumptions: recurring gravity
direction span should prevent a fixed body BA vector from maintaining its full
world projection indefinitely.

## Can gravity-direction span give a strict average contraction of held BA?

Let b_H be a fixed nonzero body-frame held accelerometer-bias vector and
`w(t)=R_true(t)' b_H` its world representation (up to the convention used in
H7). Always
`||mean w||<=||b_H||`, with equality iff w(t) has one constant direction
almost everywhere.

Gravity-direction span alone does **not** make this inequality strict
uniformly. Choose an admissible attitude motion that is a rotation about the
body axis parallel to b_H. Then w(t) is constant in world coordinates, while
the body/world gravity direction can have any prescribed nonzero span (up to
the rocking amplitude) whenever b_H is not parallel to gravity. Hence for
such histories
`Delta_g>=theta_E` but `||mean R' b_H||=||b_H||`.
Therefore
`sup_admissible ||mean R' b_H||/||b_H|| = 1`,
so no source-uniform `kappa(theta_E)<1` follows from the existing gravity-span
assumption.

Near-constant noncollinear geomagnetism does not remove this geometry by
itself: rotation about the b_H axis rotates both gravity and magnetic body
vectors and can still satisfy magnetic service. MAGNETIC SERVICE constrains
information/application, not the rotation axis relative to b_H.

Thus the hoped-for ~8% improvement cannot be obtained from Delta_g alone.
The H18 alias ceiling remains the coarse ~7.58-degree value from HO-7 unless
one uses an additional existing constraint that limits the alignment of the
held BA axis with the physical rotation axis, or exploits the full zero-action
magnetic/gyro chronology rather than only gravity span.

This is a genuine counterexample to the proposed kappa lemma, not merely a
missing estimate. Do not assert kappa(theta_E)<1 under the current MARINE
MOTION definition.

## Source covariance guard for the actual projection

The inherited BA marginal bound is `P_ba,ba <= (1/1600) I3` at every regular
operation boundary (including the pre-projection boundary). It is sharper
than the five-block nuisance comparison and is already proved in
`ou3-nuisance-upper-proof.md`. For the full, cross-coupled covariance and
`V=e'P^-1 e`, covariance Cauchy--Schwarz gives

`|e_ba| <= sqrt(lambda_max(P_ba,ba)) sqrt(V) <= sqrt(V)/40`.

With physical `B_a=.22516660498395405` and the literal projection radius .4,
projection is therefore inactive whenever, **before that projection**,

`sqrt(V) < 40(.4-B_a) = 6.993335800641838`.

In particular `sqrt(V)<=6` leaves the strictly positive distance
`.02483339501604595` from the projection boundary. The actual mean and full
covariance are then unchanged by projection, so its finite-error defect is
exactly zero. All covariance cross terms remain; Euclidean nonexpansiveness
has not been promoted to metric nonexpansiveness. This guard is uniform in
the retained source profile and independent of the missing AG upper bound.
It must be proved at every pre-projection prefix. Neither invariant radius 6,
entry into that radius, nor finite-injection/arithmetic totality follows.
The reset's injection-dependent remainder remains separately chargeable.


## Literal compatibility-graph alpha and finite-block return

The exact-kernel large-c slope can be written without choosing a soft
eigenvector. On a rank-one exact compatibility stratum use the homogeneous
literal graph

`r_W(a)=(a,0_{LIN/AW},-A_W a)`,                           (CG-1)

where `a` is a nonzero homogeneous attitude coordinate and `A_W a` is the
BA vector forced by every literal magnetic/accelerometer/S/process zero-action
compatibility equation of that SAME realized word. `A_W` therefore depends
on the actual chronological nominal-force rows, resets, held/active BA
transport and the committed coupled tuner history; it is not freely chosen.

Fix the physical metric `M>0` once and define

`s_W=r_W' M r_W`, `e_W=r_W/sqrt(s_W)`.                  (CG-2)

For a same-history successor `W_+` let `r_+` be its literal graph
generator and `n_+=M r_+/sqrt(s_+)` the corresponding unit dual functional.
On an exact-kernel word `J_W r_W=0`. In the Schur basis whose first vector
is `e_W`, positive semidefiniteness gives both the first diagonal and the
cross block of `J_W` equal to zero. Thus `j_W=0` and the Schur coefficient
in DC-3 is not arbitrary:

`ell_W=e_W' Phi_tilde_W' n_+
       = r_+' M Phi_tilde_W r_W / sqrt(s_W s_+)`.          (CG-3)

Consequently the exact large-c slope is

`alpha_W=ell_W^2
 = |r_+' M Phi_tilde_W r_W|^2/(s_W s_+)`.                 (CG-4)

If one retains the unnormalized convention of LS4, the same statement is
`alpha=ell_raw^2/s_W`; CG-4 is the invariant normalized form. On a
zero-action compatibility trajectory `Phi_tilde_W r_W=T_W r_W`, because the
conditional root map and literal deterministic homogeneous transport agree
on the data-null root. Therefore

`alpha_W
 = |r_+' M T_W r_W|^2/
   [(r_W'Mr_W)(r_+'Mr_+)]`.                               (CG-5)

This is the literal compatibility-graph transfer, not a least-eigenvector
surrogate. Every dependence of `A_W,T_W,r_+` on physical acceleration,
attitude, biases, S pseudo-updates, magnetic service, covariance-generated
gains and the coupled `tau,sigma_aw,R_S,T_S` chronology remains inside CG-5.

Hence the one-boundary theorem question is precisely

`alpha_bar_1 =
 sup_(same-history exact-kernel pairs)
 |r_+' M T_W r_W|^2/(s_W s_+) < 1 ?`.                    (CG-6)

No current lemma proves CG-6. BA decay alone cannot: the next graph can change
its attitude/BA ratio so that `T_W r_W` is collinear with `r_+`. Conversely
compactness does not prove equality; it only ensures that the maximum is
attained once the closed same-history pair class is established.

### Full-baseline block return

If `alpha_bar_1=1`, do NOT conclude instability. The correct m-word object
is the complete composed Riccati word, not the numerical composition of the
scalar functions `D_W`. Scalarization after each word discards the baseline
matrix B and its cross-covariances, which LS8--LS12 show can change the next
return.

For a same-history block
`B_m=W_(j+m-1) o ... o W_j`, compose the literal chronological factors
first, preserving the actual carried covariance, source columns, tuner state
and shared boundary states. Let

`Pi_[j,m], J_[j,m], Phi_[j,m]`

be Theorem-D's known-root covariance, root information and conditional
root-to-terminal map of that WHOLE block. With endpoint graph generators
`r_j,r_(j+m)`, define

`G_[j,m](c)=J_[j,m]+r_j r_j'/c`                         (CG-7)

in the fixed normalized physical coordinates and

`D_[j,m](c)=
 n_(j+m)' Pi_[j,m] n_(j+m)
 +n_(j+m)' Phi_[j,m] G_[j,m](c)^-1
                Phi_[j,m]' n_(j+m)`.                     (CG-8)

CG-8 is exactly LS1 applied once to the superword. It automatically includes
all intermediate quotient information and all B cross-covariance
cancellations. If the block itself has an exact endpoint kernel, its large-c
slope is

`alpha_[j,m]=
 |r_(j+m)' M T_[j,m] r_j|^2/
 [(r_j'Mr_j)(r_(j+m)'Mr_(j+m))]`,                         (CG-9)

where `T_[j,m]=T_(j+m-1)...T_j` only on the zero-action compatibility
trajectory. If any intermediate word forces positive action on the carried
mode, that mode is not in the block kernel; the block Schur information
`j_[j,m]>0` and its contribution to `D_[j,m](c)/c` tends to zero instead.

Thus a sufficient block O2 condition is: for some finite m,

`sup_(same-history exact block-kernel chains) alpha_[j,m] <= 1-delta`
for a `delta>0`.                                          (CG-10)

Together with compact full-baseline `dperp_[j,m]`, CG-10 gives a finite
block ceiling by the same LS5/DC argument. One-word unit transfer is harmless
if it cannot persist through the block.

When every constituent word is exact-kernel and the transported graph remains
exactly on each next graph, write

`T_k r_k=lambda_k r_(k+1)`.                              (CG-11)

Then CG-9 factorizes exactly:

`alpha_[j,m]=prod_(k=j)^(j+m-1) alpha_k`.                 (CG-12)

If a constituent has `alpha_k=1` but a later one has strict loss, the block
contracts. If all `alpha_k=1`, the chain is an exact persistent
compatibility execution. Therefore the finite-block alternative is equivalent
to excluding an infinite same-history unit-transfer chain, not to excluding
unit transfer on every individual word.

The existing EC compactness theorem now applies to the literal graph quantity
CG-5: if no admissible infinite recurring execution satisfies
`T_k r_k=lambda_k r_(k+1)` with unit normalized transfer at every boundary,
then diagonal compactness yields some finite `m` and `delta>0` satisfying
CG-10. Conversely an infinite equality execution defeats every such finite
block.

This is as far as the present assumptions close analytically. The coupled
physical/tuner chronology has NOT yet excluded or constructed the infinite
equality execution. Its exact equations are the already-derived compatibility
zero dynamics/BV system, now with the endpoint quantity fixed by CG-5. The
next decisive calculation is therefore to test global continuation of that
literal graph under the coupled shipping zero dynamics; further arbitrary
soft-eigenvector or separated covariance bounds cannot decide O2.


## Literal compatibility graph substituted into the coupled shipping zero dynamics

This section performs the CG substitution in the actual Live chronology.  It
does not freeze the tuner, invent an S controller, or choose an independent
soft direction.

At an accepted accelerometer epoch k let the exact homogeneous compatibility
generator be

`r_k=(a_k,0_{LIN/AW},-A_k a_k)`.                          (ZG-1)

Write `q_k=A_k a_k`.  The magnetic zero-action equations constrain `a_k`
to the transported body-field axis, while the accelerometer zero-action
equation is

`J_att,k a_k-J_ba,k q_k=0`.                              (ZG-2)

All rows in ZG-2 are the literal rows generated by the actual nominal mean,
reference, reset and covariance history.  Four-S/process zero action has
already removed an independent LIN/AW homogeneous root; it does NOT set the
base nominal AW mean to zero.

For a unit-transfer equality boundary, the homogeneous root must satisfy

`T_k r_k=lambda_k r_(k+1)`,                               (ZG-3)

with unit fixed-M normalized amplitude.  Hence, componentwise,

`a_(k+1)=lambda_k^-1 T_att,k a_k`,
`q_(k+1)=lambda_k^-1 T_ba,k q_k`.                         (ZG-4)

Substituting q=Aa gives the exact graph cocycle

`A_(k+1) T_att,k a_k = T_ba,k A_k a_k`.                  (ZG-5)

ZG-5, together with ZG-2 at every applied accelerometer/magnetic epoch, is the
literal equality manifold.  It is stronger than merely requiring a chosen
eigenvector to have unit norm.

### Shipping mean/tuner chronology on ZG

Let `z_k` denote the complete base estimator mean, covariance, tuner,
front-end, scheduler and physical state immediately before sample k.  The
shipping order is:

1. commit the schedule staged after sample k-1;
2. predict the MEKF with committed tau and sigma_aw;
3. apply a due S=0 pseudo-update with the committed actual R_S and cadence;
4. apply the accelerometer correction;
5. update the measurement-only tuner/front end from the current conditioned
   physical sample;
6. stage the new smoothed tau/sigma_aw/R_S candidate for sample k+1;
7. apply a due AW covariance sync, which is mean-neutral but changes later
   gains;
8. apply asynchronous magnetic corrections in their actual callbacks.

Thus the active schedule `theta_k=(tau_k,sigma_k,R_S,k,T_S,k)` is measurable
with respect to the preceding physical samples.  It is not an algebraic
unknown that can be chosen to satisfy ZG-5, but neither is it an extra
constraint on the current compatibility output.

Write the exact pre-accelerometer base mean after prediction and any due S
pseudo-update as

`m_k^S=S_k(theta_k,P_k) F_k(theta_k) m_(k-1)^+`,
`S_k=I-K_S,k H_S`.                                       (ZG-6)

The accepted accelerometer identity is

`u_k=r_acc,k=f_phys,k-fhat_k-bhat_a,k`,                   (ZG-7)

and the post-correction/reset map is

`m_k^+=R_k[m_k^S+K_a,k u_k]`.                            (ZG-8)

Covariance, gains and the next graph A_(k+1) follow the same actual operations.
The compatibility-maintenance equation obtained by substituting ZG-8 into the
next ZG-2 has the exact form

`h_k(z_k,q_k,u_k)=0 in R^2`.                              (ZG-9)

No independent S residual appears: it is already the state-dependent term
`-H_S m^-` in ZG-6.  No independent tuner variable appears: theta_k is the
lagged measurement-only schedule carried in z_k.

Split the physical accelerometer innovation into two transverse components
and one longitudinal component relative to the current compatibility axis,

`u_k=(u_perp,k,u_parallel,k)`.                            (ZG-10)

On every regular stratum where

`G_k=d_(u_perp) h_k`                                      (ZG-11)

is nonsingular, the implicit-function theorem gives

`u_perp,k=Psi_k(z_k,q_k,u_parallel,k)`.                   (ZG-12)

Substitution into the complete shipping map gives the exact reduced equality
dynamics

`z_(k+1)=Z_k(z_k,q_k,u_parallel,k)`,
`q_(k+1)=T_ba,k q_k/lambda_k`,                           (ZG-13)

with the tuner/front-end update and next-sample commit INCLUDED in Z_k.  This
is a one-input nonautonomous viability system on the literal compatibility
graph.

### Does the coupled tuner/S chronology force finite escape?

No.  The substitution identifies no sign-definite or divergent term.

First, the BA factor drives q toward the regular q=0 compatibility manifold;
it does not drive the trajectory toward a forbidden boundary.  At q=0,
ZG-9 is nominal field-axis compatibility and the target increment caused by
BA decay vanishes.

Second, the S pseudo-update is a bounded Kalman map inside Z_k.  It changes
the required Psi_k and the future covariance/gains but supplies no additional
independent equality beyond h_k=0.  Its innovation is endogenous and can be
bounded on a bounded compatible trajectory.

Third, the deployed tuner does not close an algebraic loop at the same
sample.  Its inputs are the measurement-only front end; its smoothed candidate
is staged after the current MEKF correction and committed before the next
sample.  Therefore its effect in ZG-13 is a bounded one-step-lagged coefficient
sequence.  The coupled law is important quantitatively--tau changes the OU
transport and T_S, sigma changes process/sync covariance, and the SpectralMSE
law changes R_S--but none of those operations adds a compatibility equation.
The safety clamps keep these coefficients in the retained compact ranges.

Fourth, AW covariance synchronization is mean-neutral.  It can alter
G_(k+1) through the next gains, so a rank-loss boundary remains possible, but
it cannot itself force the current compatible mean off ZG.

Finally, MARINE MOTION bounded velocity/displacement/potential does not create
a monotone escape supply.  After ZG-12 the free longitudinal input and the
time-varying attitude permit zero-mean periodic/quasiperiodic physical
acceleration.  Hence maintaining h=0 need not consume a nonzero DC physical
acceleration or an accumulating innovation-energy resource.

Therefore the literal coupled chronology gives

`finite escape is NOT implied by the present assumptions`.             (ZG-14)

This is a theorem about the available inequalities, not a construction of a
shipping counterexample.

### Does the substitution prove a forward-complete equality execution?

Not yet.  Local continuation follows on any strict-margin regular patch from
ZG-11--ZG-12.  A forward-complete equality execution would follow if there
were a compact positively invariant subset K of ZG on which, uniformly,

`sigma_min(G_k)>=g_0>0`,                                  (ZG-15)
all physical/retained-state and gate margins are positive, and
the recurring magnetic-service Gram has margin
`lambda_min(G_M)>=mu_M+delta_M`.                          (ZG-16)

The current theorem set proves neither ZG-15 nor invariance of those margins
along the constrained map Z_k.  In particular covariance cross terms can
drive the transverse accelerometer gain toward rank loss, and MAGNETIC
SERVICE is an assumption on the actual execution, not a proved invariant of
the compatibility-controlled continuation.

Conversely, an escape theorem would require showing every ZG trajectory
reaches one of these boundaries in finite time.  ZG-13 supplies no Lyapunov
or barrier function with such a sign.  The natural candidate |q| moves in the
wrong direction: it decays into the regular q=0 interior.

Hence the exact dichotomy after literal substitution is:

- **forced finite escape:** disproved as a consequence of the existing
  algebra/boundedness/tuner/S arguments; none supplies the required sign;
- **forward-complete compatible execution:** locally viable on regular
  strict-margin patches, globally OPEN because positive invariance of such a
  patch is not proved.

This is enough to rule out further attempts to close O2 by BA decay, S
pseudo-update accumulation, tuner coupling, bounded physical primitives, or
one-word covariance estimates.  The remaining decisive task is narrower:
construct a reachable strict-margin recurring A21 base point on ZG and prove
a compact invariant neighborhood for Z_k, OR prove a source-uniform loss of
ZG-15/ZG-16/gate margin along every constrained trajectory.  Either result
settles the infinite equality chain and therefore the finite-block O2
condition CG-10.


## Literal compatibility transport inside the full block-PSD certificate

This is the direct substitution of CG/ZG into LS3.  It is the controlling O2
object; no separately maximized d_soft or H appears.

Fix an admissible same-history m-word block
`B=[j,j+m)` and compose the literal Theorem-D factors over the WHOLE block,
with the carried physical/estimator/covariance/tuner/scheduler history.  Write

`Pi_B=Ric_B(0)`, `J_B` for the complete nuisance-shortened root information,
and `Phi_B` for the conditional root-to-terminal map.  Let the literal
endpoint compatibility graph generators be

`r_0=(a_0,0,-A_0a_0)`, `r_1=(a_1,0,-A_1a_1)`.          (BP-1)

Use one fixed physical metric M and normalize
`e_0=r_0/sqrt(r_0'Mr_0)`, `n_1=Mr_1/sqrt(r_1'Mr_1)`.
Equivalently transform root coordinates once by M^(1/2); below the Euclidean
rank-one `e_0e_0'` means that fixed metric normalization, not a wordwise
renormalization.

Define

`w_B=Phi_B' n_1`, `d_B=n_1'Pi_B n_1`,
`G_B(c)=J_B+e_0e_0'/c`.                                  (BP-2)

Then the exact block return is

`D_B(c)=d_B+w_B'G_B(c)^-1 w_B`.                          (BP-3)

All known-root process covariance, diffuse quotient uncertainty, intermediate
S/accelerometer/magnetic information, AW sync, tuner-dependent process
factors and cross-covariances are already in Pi_B,J_B,Phi_B.  In particular
Pi_B is NOT replaced by a scalar background.

The scalar ceiling is exactly the bordered PSD condition

`K_B(c):=
 [[c-d_B, w_B'],
  [w_B,   J_B+e_0e_0'/c]] >=0`.                          (BP-4)

BP-4 is LS3 for the literal superword.  Multiplying the first row/column by
sqrt(c) gives the congruent form

`Khat_B(c)=
 [[c-d_B,       sqrt(c) w_B'],
  [sqrt(c) w_B, c J_B+e_0e_0']] >=0`.                    (BP-5)

This form is useful at rank boundaries because it contains no inverse.

### Short only the true quotient, after inserting the literal graph

Choose a fixed-metric orthonormal root basis `[e_0,E_Q]`.  Write

`J_B=[[j,h'];[h,Q]]`, `w_B=(beta,g)`, `Q>0`.         (BP-6)

Here Q is the actual complete block quotient information.  Do not replace it
by an independent lower floor before taking the Schur complement.  Shorting Q
in BP-4 gives the EXACT 2x2 certificate

`S_B(c)=
 [[c-dperp,             ell],
  [ell, j+1/c]] >=0`,                                    (BP-7)

where

`dperp=d_B+g'Q^-1g`,
`j=j_B:=j-h'Q^-1h>=0`,
`ell=ell_B:=beta-h'Q^-1g`.                              (BP-8)

Thus BP-4 is equivalent to

`c>=dperp`,
`(c-dperp)(j+1/c)-ell^2>=0`.                            (BP-9)

After multiplying by c,

`j c^2+(1-ell^2-j dperp)c-dperp>=0`.                    (BP-10)

This is LS5 with unit fixed-metric root normalization, now derived after the
literal block composition.  The important point is that dperp contains the
complete background Pi_B AND the actual quotient return `g'Q^-1g`; ell
contains the correlated quotient cancellation `h'Q^-1g`.  Neither may be
bounded independently without losing the linked cancellation.

### Exact compatibility face

If the WHOLE block has a nonzero exact root compatibility mode e_0, then

`J_B e_0=0`.                                              (BP-11)

Because J_B is PSD, BP-11 forces `j=0` and `h=0` in BP-6.  Hence

`ell=beta=e_0'Phi_B'n_1
 = r_1'M Phi_B r_0/sqrt[(r_0'Mr_0)(r_1'Mr_1)]`.           (BP-12)

On the zero-action graph trajectory `Phi_B r_0=T_Br_0`, so

`alpha_B:=ell^2
 = |r_1'M T_B r_0|^2/
   [(r_0'Mr_0)(r_1'Mr_1)]`.                              (BP-13)

The complete-background certificate reduces exactly to

`S_B(c)=
 [[c-dperp_B, ell_B],
  [ell_B,     1/c]] >=0`,                                 (BP-14)

or

`c(1-alpha_B)>=dperp_B`.                                 (BP-15)

Therefore:
- if `alpha_B<1`, the exact finite ceiling is
  `c>=dperp_B/(1-alpha_B)`;
- if `alpha_B=1` and `dperp_B>0`, NO finite scalar ceiling exists for that
  block;
- if `alpha_B=1,dperp_B=0`, BP-14 is only semidefinite equality and gives no
  strict contraction.

This conclusion retains the full background.  The obstruction at alpha=1 is
not an artifact of multiplying dbar and Hbar: the positive dperp_B is the
literal known-root/quotient covariance appearing in the same Schur
certificate.

### Intermediate information and the m-word advantage

If the carried endpoint graph direction is charged anywhere inside the block,
then after complete nuisance shorting `j_B>0` unless another exact block
null mode survives.  BP-10 then has positive leading coefficient, so every
fixed block eventually satisfies BP-4 for sufficiently large c.  The positive
root is

`c_*(B)=
 [-(1-ell^2-j dperp)
  +sqrt((1-ell^2-j dperp)^2+4j dperp)]/(2j)`,             (BP-16)

with the j->0 limit given by BP-15 when alpha<1.

Thus the finite-block mechanism is sharper than multiplying boundary alphas:
a mode may have unit endpoint overlap on an early word but acquire positive
actual information later; then the superword has j_B>0 and its large-c slope
is zero.  Only a genuine exact null trajectory through the ENTIRE block lands
on BP-11--BP-15.

For a wholly exact persistent chain, ZG gives
`T_k r_k=lambda_k r_(k+1)`.  Then BP-13 factorizes into the product of the
boundary alphas under the same fixed metric.  If all are one, BP-15 fails
whenever dperp_B>0.  If some boundary is lossy while the block remains exact,
alpha_B<1 and BP-15 supplies the finite ceiling.

### Uniform same-history block theorem

For a fixed m let C_m be the compact class of admissible same-history
m-word blocks with the literal coupled tuner/physical chronology.  Define the
continuous/shorted BP quantities on each closed event stratum.  A uniform
block ceiling exists if and only if the BP-10 positive roots are uniformly
bounded.  A sufficient and, on the exact-kernel face, necessary condition is

`sup_(B in C_m: j_B=0) alpha_B <1`,                      (BP-17)

together with the already required compact finite `dperp_B` and positive
quotient shorting on the retained coordinates.  Then choose

`c_m >= sup_(B in C_m) c_*(B)<infinity`                  (BP-18)

and BP-4 holds for every block.  Strict inequality can be retained by choosing
c_m above the attained supremum and preserving the existing nonlinear/
arithmetic margins.

If BP-17 fails for every finite m because there is an infinite exact
unit-transfer compatibility execution, the scalar block ceiling architecture
cannot close: every prefix lies on BP-15 with alpha=1 and positive background.
If no such infinite execution exists, the EC compactness theorem supplies
some finite m and delta>0 on the exact face; continuity of BP-10 then gives a
finite uniform c_m.  This is the precise bridge from the compatibility
continuation problem to the full-baseline linked Riccati certificate.

No theorem flag is promoted here: the existence/nonexistence of the infinite
exact equality execution remains OPEN.  What is closed is the algebraic
question of how its answer enters O2: through BP-4/BP-10, with the complete
background covariance retained.


## Can q=0 compatibility acquire positive full-block information?

This calculation addresses the remaining possibility after BP: perhaps the
endpoint graph transfer stays unit but the complete block nevertheless charges
the carried mode, so `j_B>0`.

Use the exact fixed-factor block model after all literal operations have been
composed:

`y=O_s x+A xi+O_f x_f`,                                   (QI-1)

where xi stacks every fresh process, measurement and sync source ONCE and x_f
is the nuisance root.  With the actual joint source covariance whitened into
the factor columns, nuisance shorting gives

`x'J_B x =
 min_(x_f,xi) ||O_s x+O_f x_f+A xi||^2+||xi||^2`           (QI-2)

in the equivalent minimum-action representation (or its exact correlated
whitened form).  The essential fact is positivity: `x'J_Bx=0` iff there is
a literal homogeneous trajectory with ZERO total source/measurement action.
No cancellation between positive source actions can make a nonzero action
zero.

Let `e_0` be the normalized q=0 compatibility root.  The shorted scalar in
BP-8 obeys

`j_B=0 <=> e_0'J_Be_0=0`                                 (QI-3)

because J_B is PSD and the quotient Schur complement is its minimum over the
quotient coordinates.  Thus j_B is exactly the minimum complete-block action
of the carried graph mode after all allowed nuisance-root mimics, not a
separate observability constant.

### Zero-action classification on q=0

Set q=0 in the literal graph.  If `j_B=0`, every nonnegative component of
QI-2 must vanish.

1. **Fresh process/sync factors vanish.**  AW, LIN, AG/BA fresh process factors
   and AW covariance-sync factors are zero in the homogeneous auxiliary
   trajectory.  Hence no intermediate nuisance can be reselected from word to
   word; it is the deterministic image of the block root.

2. **S rows vanish.**  Every applied S=0 row has zero homogeneous residual.
   With zero fresh LIN/AW source action, the four-S injectivity result removes
   the independent homogeneous `(v,p,S,a_w)` root on regular blocks.  This
   statement concerns the error trajectory; the BASE nominal AW mean remains
   free to follow the coupled tuner/physical zero dynamics.

3. **Magnetic rows vanish.**  Every applied informative magnetic row has zero
   homogeneous residual.  MAGNETIC SERVICE plus the deterministic AG transport
   leaves the transported field-axis attitude class (and the already-qualified
   gyro-bias restrictions).

4. **Accelerometer rows vanish.**  With q=0 and no homogeneous LIN/AW/BA mimic,
   the homogeneous accelerometer condition is exactly the field-axis
   compatibility equation generated by the ACTUAL base nominal specific-force
   history.  In the notation of ZG this is `h_k(z_k,0,u_k)=0` at every
   applied epoch.

Therefore

`j_B=0`
`<=>`
`the q=0 carried graph extends as a zero-action compatibility trajectory
 through every literal event of B`,                                    (QI-4)

subject to the already stated four-S/nullity qualifications.  The reverse
direction is immediate: such a literal zero-action trajectory is an admissible
competitor in QI-2 with zero cost.

This equivalence is stronger than the earlier endpoint alpha statement.
Endpoint unit transfer alone does NOT imply j_B=0; any intermediate violation
of S/magnetic/accelerometer compatibility gives positive block action.  But an
exact ZG trajectory does imply j_B=0 for every finite prefix.

### Does the coupled q=0 zero dynamics necessarily make j_B positive?

No theorem in the present assumptions does so.  Substitution of ZG-12 into
the base shipping recursion shows why.  On every regular strict-margin patch,
the two transverse physical accelerometer components solve the two
compatibility equations:

`u_perp,k=Psi_k(z_k,0,u_parallel,k)`.                     (QI-5)

The S pseudo-update is already inside z_k and Psi_k.  The lagged
`tau,sigma_aw,R_S,T_S` tuner schedule changes the coefficients and future
covariance/gains but adds no homogeneous action when the corresponding
homogeneous S/process residuals are zero.  Magnetic service constrains the
homogeneous attitude mode but is compatible with the field-axis line.  Hence
none of these terms creates an unavoidable positive summand in QI-2 while the
constrained base trajectory remains on a regular ZG patch.

In particular there is no valid implication

`endpoint alpha=1 => j_B>0 after at most m words`          (QI-6)

from the current assumptions.  Proving QI-6 would be exactly an escape theorem
for ZG in different notation.

Conversely, local IFT viability is not a proof that j_B=0 can persist forever:
the constrained base trajectory may eventually hit transverse-rank loss,
a physical/gate boundary, or the magnetic-service boundary.  At that first
event the exact zero-action continuation can fail and every sufficiently long
superword containing the failure has `j_B>0` for the carried mode.

### Compactness equivalence for finite-block rescue

Let ZG_infty be the inverse-limit class of admissible recurring q=0
same-history executions.  Assume the established compact event-stratum and
four-S/nullity qualifications.  Then the following are equivalent:

(A) there is no forward-complete execution in ZG_infty satisfying exact
    compatibility at every literal event;

(B) there exist finite m and epsilon_J>0 such that every admissible m-word
    q=0 carried graph mode has
    `j_B>=epsilon_J` OR leaves the exact block-kernel face with
    `alpha_B<=1-delta` for some uniform delta>0.           (QI-7)

Proof.  If (B) fails for every m, choose longer and longer blocks whose
shorted action tends to zero while endpoint loss tends to zero.  Compactness,
lower semicontinuity of the nonnegative joint action and diagonal extraction
produce an infinite literal zero-action ZG execution.  Conversely an infinite
exact ZG execution has j_B=0 and alpha_B=1 on every finite prefix, contradicting
(B).

A slightly stronger pure-information statement,
`j_B>=epsilon_J>0` on every m-word block, is NOT equivalent and need not
hold: a finite block can remain an exact kernel but have alpha_B<1, which is
already sufficient through BP-15.  The correct rescue alternative is the
union in QI-7: positive block information OR strict exact-kernel transfer loss.

### Consequence for the proof strategy

The full-block information side does not independently eliminate the equality
execution.  It proves that the information and compatibility formulations are
the SAME obstruction:

`infinite exact q=0 ZG execution`
`<=>`
`for every finite prefix: j_B=0 and alpha_B=1`.           (QI-8)

Thus searching for a generic positive lower bound on j_B while allowing exact
ZG compatibility would be circular.  The next useful theorem must operate on
the BASE constrained dynamics, not on another homogeneous information bound:
either prove every q=0 ZG trajectory reaches a rank/service/gate/physical
boundary in finite time, or construct one compact forward-invariant
strict-margin ZG execution.  Once that is decided, BP-10 converts the result
directly into (or rules out) finite-block O2.


## Quantitative transverse authority on the q=0 compatibility manifold

This calculation differentiates the literal next compatibility output with
respect to the CURRENT accepted accelerometer innovation.  It keeps the
shipping covariance/gain chronology; G is not replaced by an AW-gain proxy.

Let x_k^S be the estimator mean immediately before accepted accelerometer
correction k, after prediction and any due S pseudo-update.  Let C_k be the
literal 3-row accelerometer Jacobian used by shipping,

`C_k=[J_att,k, 0_bg, 0_vps, R_wb,k, J_ba,k]`              (GA-1)

with the optional lever-arm gyro-bias block inserted when enabled and with the
BA block omitted from the gain when BA updates are frozen.  Shipping computes

`S_k=C_k P_k C_k'+R_acc,k`,
`K_k=P_k C_k' S_k^-1`.                                   (GA-2)

The mean increment is K_k u_k, followed by the literal quaternion injection,
bias projection (on a strict interior patch its derivative is identity), and
error-state reset.  Denote the derivative of that complete post-correction
mean/reset map by R_k.  Compose from there to the next compatibility evaluation
all literal predictions, magnetic corrections, and any due S pseudo-update;
call this derivative F_(k+1,k).  It uses the committed lagged
`tau,sigma_aw,R_S,T_S` schedule and the covariance-generated gains on that
same base history.

At q=0 the next homogeneous compatibility output is the two-dimensional
projection of the next nominal-force/attitude row transverse to the
transported magnetic axis.  Let C_c,k+1 be its derivative with respect to the
base mean and let E_perp,k inject the two selected transverse components of
the physical accelerometer innovation.  Then, away from gate/projection
boundaries,

`G_k=d_(u_perp) h_k
 =C_c,k+1 F_(k+1,k) R_k K_k E_perp,k`                     (GA-3)

plus the direct physical-input term if h_(k+1) is defined at the same sample.
For the one-step-ahead convention used in ZG there is no direct term.  In
world field-axis notation one may write
`C_c=P_bperp C_m`, giving the same formula.

Substitute GA-2:

`G_k=N_k S_k^-1 E_perp,k`,                                (GA-4)
`N_k:=C_c,k+1 F_(k+1,k) R_k P_k C_k'`.                   (GA-5)

Since S_k is SPD on every accepted update, rank(G_k) is exactly the rank of
the corresponding two-column restriction of N_k.  Positive R_acc cannot
supply transverse authority that is absent from N_k.

### Exact singular-value comparison

Let `S_perp,k=E_perp,k' S_k E_perp,k` only when E_perp selects an invariant
measurement plane.  In the general case keep the rectangular right factor.
For any 2-vector z,

`|G_k z| >= sigma_min(N_k S_k^-1 E_perp,k)|z|`.           (GA-6)

Using singular-value products gives the valid coarse lower implication

`sigma_min(G_k)
 >= sigma_min(N_k|Range(S_k^-1 E_perp,k))
    sigma_min(S_k^-1 E_perp,k)`.                          (GA-7)

On the retained compact covariance/noise class,
`sigma_min(S_k^-1 E_perp,k)>=1/lambda_max(S_k)>0`.
Therefore a uniform authority floor is EQUIVALENT, up to known finite
conditioning, to a positive source-uniform floor for the transported
cross-covariance numerator N_k on the actual transverse innovation plane.

The coupled tuner law supplies compact upper/lower bounds for the conditioning
factor through process covariance, R_S and cadence, but it does not by itself
give a lower singular bound for N_k.

### Why covariance positivity does not force N_k to be nonsingular

The numerator is a cross covariance between the current accelerometer
measurement and the NEXT compatibility output after the intervening corrected
chronology:

`N_k=Cov(h_(k+1), y_acc,k | past)`                        (GA-8)

in the linearized joint Gaussian model (with the literal reset/transport).
A positive-definite state covariance P_k and positive R_acc guarantee
S_k>0, but a cross covariance may vanish.

This is not merely a loose-bound issue.  Partition the current state into the
two-dimensional compatibility-output sector c and the remaining state r.
Then the relevant numerator has the schematic exact form

`N_perp=A P_cc H_c'
        +A P_cr H_r'
        +B P_rc H_c'
        +B P_rr H_r'`.                                    (GA-9)

The off-diagonal covariance blocks are signed.  Joseph updates, S
pseudo-updates, magnetic updates and predictions preserve PSD of the WHOLE
P but do not preserve the sign or a lower singular value of this particular
cross block.  AW covariance sync in the deployed path adds a PSD increment to
the AW marginal while preserving existing cross-covariances; this changes
relative correlations but again supplies no sign constraint on GA-9.

A two-state SPD witness already shows the algebra: with
`P=[[1,rho],[rho,1]]`, measurement row H=[1,0], and next output row
L=[-rho,1], one has `L P H'=0` for every |rho|<1 although P>0 and
S=1+R>0.  This witness is algebraic, not claimed shipping reachable.  It
proves that covariance positivity/noise floors alone cannot establish the
desired authority floor.

### What the actual tuner/S chronology does and does not guarantee

The committed `tau` bounds keep OU prediction coefficients finite and away
from their singular limits on each positive-dt step.  Positive sigma_aw and
the pending AW floor provide process covariance in the AW sector.  The
SpectralMSE R_S law and bounded T_S give recurring finite-noise S corrections.
Together these are valuable for compactness and for bounding P and S_k.

But G_k depends on the ORIENTATION of the full covariance through P_k C_k'
and on its subsequent transport through F R.  None of the coupled scalar
laws fixes that orientation.  S and magnetic Joseph corrections can rotate
the relevant cross-covariance; accelerometer corrections can do the same.
Therefore the coupled law does not imply

`inf_(q=0 compatible histories) sigma_min(G_k)>0`.         (GA-10)

Conversely it also does not imply inevitable rank loss.  Full rank is an open
condition: if one reachable q=0 A21 point has `det G_k !=0`, then a
neighborhood of that point has a positive local floor.  Constructor/diagonal
covariance examples establish algebraic rank-two authority, but they are not
yet certified reachable recurring A21 q=0 roots.

### Indefinite persistence versus inevitable loss

Define the regular authority set

`R_g={z in ZG(q=0): sigma_min(G(z))>=g}`.                 (GA-11)

For any g>0, R_g is closed inside a fixed event/gate stratum; the strict set
sigma_min(G)>g is open.  The shipping constrained map Z sends a regular point
to its next q=0 point after solving the transverse innovation.

The current equations prove neither

`exists g>0, compact K subset R_g with Z(K) subset K`      (GA-12)

nor

`every q=0 constrained trajectory reaches det G=0 in
 finite time`.                                             (GA-13)

Thus the answer to the proposed question is precise:

**the coupled covariance/tuner/S chronology CAN remain uniformly nonsingular
only if a compact invariant regular set GA-12 exists; nothing in the current
scalar tuning laws forbids it.  It is NOT forced to remain nonsingular, and
no current invariant forces eventual singularity.**

This means rank loss cannot presently be used as the automatic mechanism
that produces j_B>0.  If a constrained trajectory reaches det G=0, exact
compatibility may still continue through that point using a different input
chart; singularity of this particular 2x2 chart is not itself positive
homogeneous information.  To infer j_B>0 one must prove that the FULL
compatibility equation has no admissible continuation there, not merely that
one transverse parameterization fails.

### Stronger chart-invariant authority criterion

Let D_k be the full derivative of the two compatibility equations with
respect to the THREE physical accelerometer components:

`D_k=C_c,k+1 F_(k+1,k) R_k K_k`.                         (GA-14)

Local compatibility continuation requires `rank D_k=2`; a choice of
E_perp is only a coordinate chart.  Therefore the intrinsic quantity is

`g_full(z)=sigma_2(D_k)`,                                 (GA-15)

the second singular value.  There exists some transverse 2-plane with a
well-conditioned G iff g_full>0.  If one selected G becomes singular while
g_full>0, change charts; no compatibility break has occurred.

Using K=P C'S^-1,

`D_k=C_c F R P C' S^-1`.                                 (GA-16)

As before, S^-1 is uniformly conditioned on the retained class, so the
decisive quantity is the rank-two cross-covariance operator

`C_c F R P C'`.                                          (GA-17)

A finite-block information gain can be forced from authority loss only after
proving `rank D_k<2` makes the affine compatibility equation unsolvable for
the actual drift term.  Rank loss alone may instead leave a rank-one or
rank-zero compatible solution if the drift lies in Range(D_k).

Therefore the next quantitative target should not be a lower bound on one
chosen G.  It is the chart-invariant constrained solvability margin

`eta_k=dist(-h_k(z,0,0), Range(D_k))` when rank D_k<2,     (GA-18)

together with `sigma_2(D_k)` on the regular set.  A source-uniform theorem
that every infinite q=0 constrained execution either stays in a compact
`sigma_2(D)>=g>0` set or incurs `eta>=eta_0>0` at a rank-loss event would
decide continuation versus positive block information without coordinate
artifacts.

No such theorem is currently proved.  This calculation rules out the simpler
hope that positive covariance, R_acc, S recurrence, or the coupled tuner law
alone enforce or destroy transverse authority.


## Shipping-closed self-consistency equation for the alleged q=0 trajectory

The preceding ZG/GA formulations deliberately exposed local compatibility
authority, but they still permit a misleading reading: choose an innovation
u to keep h=0 and ask later whether some physical motion realizes u.  Shipping
does not have that causal freedom.  The accelerometer residual, the private
measurement-only front end, the tuner schedule, covariance, Kalman gain and
next residual are generated by ONE physical IMU history.

This section therefore removes u as an independent control.

### One literal sample map

Let the exogenous physical sample be

`p_k=(R_k^phys, omega_k^phys, a_k^phys, b_k^phys, temperature_k, ...)`

with calibrated accelerometer sample `y_k=Y(p_k)`.  Let `z_k` contain the
complete shipping state immediately before updateCore_: MEKF mean/covariance,
private Mahony/front-end state, wave-period/frequency/variance state, tuner
EMAs and targets, pending one-sample tune commit, S/AW-sync clocks, magnetic
reference/service state, and gate state.

Shipping executes one deterministic map

`z_(k+1)=S_k(z_k,p_k)`.                                   (SC-1)

Its relevant factorization is causal:

`theta_k = Commit(z_k.pending)`,                          (SC-2)
`ybar_k = conditionAccel(y_k; z_k)`,
`front_(k+)=Front(front_k,omega_k,ybar_k)`,               (SC-3)
`(x_k^-,P_k^-)=Predict_theta_k(x_k,P_k)`,
`(x_k^a,P_k^a)=AccUpdate(ybar_k;x_k^-,P_k^-)`,            (SC-4)
`tuner_(k+)=Tune(front_(k+),ybar_k;tuner_k)`,
`pending_(k+)=Stage(tuner_(k+))`,                         (SC-5)
followed by due AW covariance sync and asynchronous magnetic callbacks in
their literal order.  A due S=0 pseudo-update is inside Predict_theta_k with
the committed actual R_S and T_S.  The exact source order in the wrapper is
commit -> condition -> private level/front end -> MEKF prediction/acc update
-> tuner -> AW sync; the period estimator used by tuning is lagged as in the
shipping code.

The committed parameter tuple is therefore a functional of the PRIOR physical
history,

`theta_k=(tau_k,sigma_aw,k,R_S,k,T_S,k)
          =T_k[p_0,...,p_(k-1);z_0]`.                     (SC-6)

For the deployed SpectralMSE branch this functional includes the
measurement-only period-scaled acceleration band, variance horizon, operating
point
`tau_target=c_tau T_z/2`,
`sigma_target=c_sigma sigma_a,B`,
the tau/sigma EMA, and

`R_S,target=C_J q_eff^(1/14) sigma_a,B^(6/7)
              tau^(24/7)/sqrt(T_S)`,                      (SC-7)

with `T_S` the actual tau-scaled/clamped pseudo cadence, R_S EMA/clamps, and
the one-sample delayed commit.  Thus tau,sigma_aw,R_S,T_S are not independent
unknowns in the compatibility equations.

### Impose q=0 compatibility BEFORE solving for the physical input

On the alleged late pathological branch, the homogeneous compatibility line
has q=0 and the nominal predicted specific force is parallel to the
transported magnetic field.  In the fixed-world form already derived,

`P_bperp a_hat_w,k = P_bperp g + delta_field/lever,k`.     (SC-8)

Call the right-hand required nominal AW value `a_hat_w,k^*(z_k,p_k)`;
its transverse component is gravity scale.  The literal pre-update prediction
and any due S correction produce `a_hat_w,k^-`.  Shipping accelerometer
correction gives

`a_hat_w,k^+
 =a_hat_w,k^- + K_aw,k(z_k,theta_k,ybar_k) r_k`,           (SC-9)

where

`r_k=ybar_k-h_acc(x_k^-,P_k^-,temperature_k)`.            (SC-10)

Attitude, BA and lever terms are corrected by OTHER rows of this SAME K_k and
the quaternion/reset map, so the exact compatibility equation is more
generally

`C_comp,k S_k(z_k,p_k)=0`.                                (SC-11)

When the literal correction map is locally invertible on the two required
compatibility components, SC-11 determines the residual required by
compatibility as a FUNCTION of the pre-sample shipping state and the remaining
physical sample components,

`r_k=r_req(z_k,p_k;theta_k)`,                              (SC-12)

not as a free control.  Substituting the measurement identity SC-10 gives the
required physical accelerometer sample

`ybar_k =
 h_acc(x_k^-,temperature_k)+r_req(z_k,p_k;theta_k)`.       (SC-13)

After undoing the deterministic conditioning/calibration map on its accepted
branch, SC-13 is an equation for the physical acceleration:

`a_k^phys =
 A_req[z_k,R_k^phys,omega_k^phys,temperature_k;theta_k]`.  (SC-14)

### Close the tuner/covariance loop

But theta_k and the covariance/gain appearing in A_req are themselves
generated by the SAME physical history.  Substitute SC-6 into SC-14 and
advance the complete shipping state with SC-1:

`a_k^phys =
 A_req[z_k,R_k^phys,omega_k^phys,T_k[p_<k;z_0]]`,          (SC-15)
`z_(k+1)=S_k(z_k,p_k)`.                                   (SC-16)

Equivalently, over an m-word physical history p define

`T[p]` = literal front-end/tuner schedule functional,
`C[p]` = literal Riccati/Joseph/AW-sync covariance functional driven by T[p],
`R[p;T[p],C[p]]` = residual sequence required by q=0 compatibility,
and `A[p;...]` = physical acceleration reconstructed from that residual.

The pathological trajectory exists only if

`a_phys = F_ship[a_phys]
 := A[R[a_phys;T[a_phys],C[a_phys]],
      T[a_phys],C[a_phys]]`                                (SC-17)

together with the physical attitude/gyro/magnetic history and every MARINE
MOTION, IMU BIAS, MAGNETIC SERVICE, gate and retained-state condition.

SC-17 is the requested shipping self-consistency equation.  It is a
history-dependent delayed nonlinear fixed point, not the pointwise
two-equation/two-control viability problem.

### Immediate analytical consequences

1. **Independent extrema are forbidden.**  Any proof step choosing a
   compatibility-favorable residual, a different tuner-favorable physical
   spectrum, and an independently favorable covariance is invalid.  They must
   all be images of the same p under SC-1.

2. **The gravity-scale nominal AW requirement is not a physical acceleration
   requirement.**  SC-8 forces a persistent nominal AW transverse component,
   while bounded physical velocity forces zero long-time mean physical
   acceleration.  Therefore any fixed point must continually regenerate the
   gravity-scale nominal AW through the gain-weighted residual/S/magnetic
   correction supply while its physical acceleration remains AC.  This is the
   exact coupled balance to test.

3. **The tuner sees physical/front-end acceleration, not r_req.**  In
   particular the variance channel uses the private levelled acceleration,
   period-scaled band and lagged wave-period estimate.  A large DC residual
   caused by nominal-force mismatch need not appear as tuner wave variance in
   the same way; conversely an oscillatory physical acceleration that makes
   the tuner choose a given sigma/tau also changes the residual in SC-10.
   These effects cannot be separated.

4. **Covariance is deterministic once the physical/event history is fixed.**
   Given z_0 and p, Joseph updates, S cadence/R_S, OU Q(tau,sigma), magnetic
   updates and AW sync determine P and K.  The gain needed in SC-12 cannot be
   selected independently.

5. **Commit delay matters.**  The sample y_k that must satisfy SC-13 cannot
   alter theta_k used by its own correction.  It only changes the candidate
   committed at k+1.  Hence a putative periodic/recurrent fixed point must
   close the augmented delayed state, including tuner/front-end memory and
   pending commit, not merely the MEKF mean.

### Eliminate the residual exactly

On an accepted branch, SC-10 lets us remove r entirely from the fixed-point
unknowns.  Define the compatibility residual map

`H_k(z_k,p_k):=C_comp,k S_k(z_k,p_k)`.                    (SC-18)

The exact pathological history is simply a physical history satisfying

`H_k(z_k,p_k)=0` for every required accelerometer/magnetic epoch,
`z_(k+1)=S_k(z_k,p_k)`,                                   (SC-19)

with theta/P/K generated internally by S.  This formulation is tautological
but important: it prevents the proof from enlarging the reachable set by
promoting innovations, gains, covariance or tuner parameters to controls.

For analysis, split p into the actual physical acceleration a and the
remaining physical history w.  The local derivative is the TOTAL shipping
derivative

`D_ship,k=d_a H_k
 =partial_a H_k
  +partial_z H_k * d_a z_k
  +partial_theta H_k * d_a T[p_<k]
  +partial_P H_k * d_a C[p_<k]`.                          (SC-20)

At the current sample, causality makes the last three history terms depend
only on prior acceleration samples; over a word they are essential.  The
earlier GA operator retained only the instantaneous correction authority with
the base history frozen.  SC-20 is the correct derivative for reachability of
the shipping fixed point.

### What is proved now, and what is not

SC-17/SC-19 does NOT yet prove that the pathological fixed point is impossible.
It proves that the generic zero-dynamics viability calculation is only a
relaxation and cannot establish reachability of the pathology.

Nor can the existing local rank-two gain result establish existence: it proves
solvability after freezing a base history, whereas SC-17 requires that the
resulting physical samples reproduce that same base history's tuner and
covariance sequence.

The next falsifiable analytical calculation is now well posed.  Assume an
asymptotically periodic/recurrent physical history of period L on q=0.  Lift
SC-1 over L to the augmented state including tuner/front-end/covariance memory:

`z_L=Phi_L(z_0,p_[0,L))`,
`H_[0,L)(z_0,p)=0`.                                      (SC-21)

A recurrent pathology requires a fixed point (or compact recurrent orbit) of
this SAME map, with the physical acceleration satisfying its zero-mean and
moment constraints.  One should next eliminate the MEKF linear mean states
from SC-21 using the exact OU/S transition for the tuner sequence generated by
p, leaving a reduced fixed-point equation in the physical waveform and the
front-end/tuner/covariance orbit.  If the resulting equations are inconsistent,
the alleged pathology is unreachable and finite-block O2 follows by
compactness.  If they have one strict-margin admissible solution, that is a
genuine shipping-reachable counterexample to this proof architecture.

Until SC-21 is solved or excluded, no generic viability trajectory should be
called an admissible obstruction.


## Candidate replacement proof: complete-word dissipativity / LaSalle feasibility test

This section is deliberately parallel to O1/O2. It tests whether the exact
source-faithful complete-word factorization already supplies a useful
semidefinite storage/dissipation identity before any scalar kernel ceiling,
diameter, generalized eigenvalue or independent extremum is introduced.

### Exact minimum-action storage

Freeze one literal recurring linear A21 service word W, including the actual
same-history tuner/covariance/event chronology. Stack the root x, every fresh
whitened process/sync/noise factor s exactly once, all applied measurement
residuals y, and terminal error x_+. The existing complete-word factorization
has the form

y = O x + A s,
x_+ = T x + B s.                                           (DL-1)

The whitening uses the actual joint source covariance, so correlations are
retained. Split x=(x_s,x_f) only when eliminating a genuine nuisance root.

For a prescribed root x define the complete zero-output action

D_W(x) :=
 min_(x_f,s) { ||s||^2 + ||R^(-1/2)(O_s x+O_f x_f+A s)||^2 },  (DL-2)

where the second norm is shorthand for the exact whitened applied
accelerometer, magnetic and S rows; equivalently absorb R^(-1/2) into O,A.
This is the same fixed-factor variational object used by the complete-word
information proof. Therefore

D_W(x)=x' J_W x >=0.                                       (DL-3)

No covariance ceiling or O1/O2 scalar bound is required for DL-3.

The dynamic-programming form is more revealing. Let V_k(e) be the minimum
remaining whitened source/measurement action from event phase k to the end of
the frozen word, conditional on current homogeneous error e. Then every
literal operation satisfies the Bellman equality

V_k(e_k)=min_(fresh source)
 { d_k(e_k,source)+V_(k+1)(e_(k+1)) },                     (DL-4)

with d_k>=0 the exact local process or innovation action. Along the minimizing
trajectory,

V_k(e_k)-V_(k+1)(e_(k+1))=d_k.                             (DL-5)

Summing the actual allowed event path gives

V_root(x)-V_terminal(x_+)=D_W(x)
 =D_proc+D_sync+D_S+D_acc+D_mag+D_terminal >=0.            (DL-6)

DL-6 is the requested path-complete dissipativity identity. The phase
storages V_k are quadratic value functions (Schur complements of the same
joint Gaussian action); they need not decrease under every physical operation
when viewed in one fixed Euclidean metric. Their Bellman differences are
nonnegative on the allowed event graph.

For a disturbance/noise input w that is not minimized as an internal source,
completion of squares gives the supply form

V_+-V_- <= -D_W(e)+w' Q_W w                               (DL-7)

after augmenting the value function with the exact cross term, or the
equivalent joint quadratic supply matrix before Young relaxation. The
homogeneous feasibility test needs only DL-6.

### Equality set: no scalar ceiling is needed

Because DL-6 is a sum/minimum of nonnegative whitened actions,

D_W(x)=0                                                   (DL-8)

if and only if one literal homogeneous trajectory makes every constituent
action zero simultaneously. Hence equality forces:

1. every fresh process and AW-sync source factor is zero;
2. every applied S pseudo-measurement homogeneous residual is zero;
3. every accepted accelerometer homogeneous residual is zero;
4. every applied magnetic homogeneous residual is zero;
5. terminal/minimum-action nuisance loss is zero.

These are not independently chosen zeroes: the SAME deterministic homogeneous
trajectory must satisfy all five through the literal chronological maps.

Existing proof lemmas then apply without quantitative scalarization.
Four-S injectivity removes an independent homogeneous LIN/AW root on a regular
service word. Magnetic zero action plus service restricts AG to the transported
field-axis class. Zero fresh BA action makes BA follow its literal homogeneous
transport. Accelerometer zero action intersects those classes in the physical
tilt/BA compatibility graph. Thus, subject to the already recorded four-S and
multi-epoch intersection qualifications,

Ker D_W = Ker J_W subseteq C_W,                             (DL-9)

where C_W is the literal compatibility graph. On strata where the physical
compatibility line is known to exist exactly, equality gives that line rather
than strict one-word dissipation.

This is a GOOD outcome for the proposed architecture: the feasibility test
does not need to prove one-word strictness.

### Why this is not merely O1 in new notation

O1/O2 next ask for a quantitative lower eigenvalue away from the kernel and a
separate covariance ceiling along it. DL instead keeps only the exact
semidefinite action and follows its equality trajectory through successive
same-history service words.

Let e_(j+1)=F_Wj^(0)e_j denote the deterministic homogeneous transport selected
by D_Wj(e_j)=0. Define the m-block accumulated dissipation

D_[j,m](e_j)=sum_(r=0)^(m-1) D_W(j+r)(e_(j+r)).             (DL-10)

Then D_[j,m]=0 iff the SAME error trajectory lies in the literal zero-action
set of every constituent word. No intermediate compatibility direction,
covariance, tuner tuple or nuisance mimic may be reselected independently.

The central LaSalle question is therefore

largest invariant subset of {D_W=0} under literal shipping chronology = ?  (DL-11)

This is exactly the shipping-closed compatibility fixed-point problem derived
in SC/PE. If the alleged compatibility history is not reachable by the actual
coupled physical/front-end/tuner/covariance equations, it is not in the
invariant equality set even though each relaxed word separately has a
compatibility graph.

### Compactness gives strict finite-window dissipation

Assume a compact retained recurring A21 execution class, continuity/lower
semicontinuity of the fixed-factor action across its finite event strata, and
uniform equivalence m0|e|^2 <= V_phase(e) <= m1|e|^2. Suppose the only
forward-complete same-history trajectory satisfying D_Wj(e_j)=0 for every j
is the desired zero error (or explicitly accepted gauge).

If no finite-window strictness existed, for every n there would be a normalized
initial error and n-word same-history execution with accumulated dissipation
tending to zero. Event compactness and the fixed-factor lower-semicontinuity
argument already used elsewhere in the proof give a diagonal limiting
infinite trajectory with D_Wj=0 for every j, contradicting DL-11.

Therefore there exist finite m and epsilon>0 such that

sum_(r=0)^(m-1) D_W(j+r)(e_(j+r))
 >= epsilon V_j(e_j).                                      (DL-12)

Telescoping DL-6 gives

V_(j+m)(e_(j+m)) <= (1-epsilon) V_j(e_j).                  (DL-13)

Thus strict contraction follows from invariant-set exclusion plus compactness,
without a kernel covariance ceiling.

For disturbances, DL-7 and the same m-block sum give regional ISS/practical
stability once the nonlinear/reset remainder is absorbed by the existing
small-radius estimates.

### Does the experiment pass?

At the analytical level: YES, conditionally.

- A nonnegative complete-word quadratic action already exists in the current
  source-faithful factorization.
- It admits a path-complete Bellman/storage interpretation.
- Its equality conditions are exactly the simultaneous zero process/S/acc/mag
  actions already characterized by the proof.
- No scalar kernel ceiling is needed to state nonexpansiveness or the equality
  set.
- The remaining strictness problem is precisely the same-history invariant
  zero-dissipation problem, where the new shipping-closed tuner/covariance
  self-consistency analysis belongs.

What is NOT yet proved is the decisive DL-11 invariant-set exclusion. Nor has
a single globally uniform phase metric/value function been constructively
enclosed over every shipping-generated tuner/covariance history. Those are the
two obligations to test before replacing O1/O2 as the main theorem path.

The experiment therefore justifies continuing this architecture, but not yet
deleting O1/O2.

### Next decisive calculation

Use SC-19/PE on an equality trajectory, but now equality supplies MORE than
q=0 compatibility: every homogeneous fresh source, S, accelerometer and
magnetic action is exactly zero. Substitute those equality conditions into
the shipping-closed physical/front-end/tuner/covariance recursion and ask
whether a nonzero forward-complete same-history trajectory exists.

If none exists, DL-12 follows by compactness and the scalar O1/O2 ceiling path
can be retired. If one exists, it is a genuine zero-dissipation obstruction
to this stronger architecture and immediately identifies the theorem/assumption
gap.


## Periodic base field-alignment Fredholm test for the LaSalle obstruction

This calculation tests the single surviving DL equality candidate on the BASE
shipping system. It imposes

P_B a_w,k = P_B g                                             (PF-1)

at every required accelerometer epoch, where P_B projects onto the fixed-world
plane transverse to the qualified geomagnetic reference (with the declared
field/lever defect tube suppressed here for readability). The physical
waveform is the only external degree of freedom. Tuner parameters, covariance,
gains and event chronology are generated by that same waveform.

Fix an L-periodic physical waveform p and suppose its complete private-front-
end/tuner/covariance/event orbit rho is also L-periodic. Conditional on this
same-history rho, collect the affine MEKF mean coordinates in x. The literal
shipping recursion is

x_(k+1)=A_k(p,rho)x_k+b_k(p,rho).                            (PF-2)

Lift one period:

x_L=A_L x_0+b_L.                                            (PF-3)

Stack every PF-1 row as C_k x_k=d_g. After substitution of the lifted
intermediate states this gives

H_B x_0 = q_B.                                              (PF-4)

Hence periodic field-aligned base means exist iff the finite Fredholm system

M_PF x_0=q_PF,                                              (PF-5)
M_PF=[I-A_L; H_B],  q_PF=[b_L; q_B]

is consistent. Equivalently every left-null vector l of M_PF must satisfy
l' q_PF=0. This is the exact mean-state arithmetic test; no innovation or gain
is promoted to a control.

If rho(A_L)<1, the mean is unique and PF-5 reduces to

F_PF(p,rho):=H_B(I-A_L)^(-1)b_L-q_B=0.                      (PF-6)

If A_L has unit Floquet modes, PF-5 is retained without regularization.

### Does the literal AW/LIN/S structure make PF-5 inconsistent?

No. The existing frozen S-to-S calculation is the one-direction Schur
reduction of PF-5. For a transverse direction, with lifted mean map
x_+=A_S x+B_S u and AW selector e_a, compatibility e_a'x=g_perp gives

G_DC=e_a'(I-A_S)^(-1 B_S,                                  (PF-7)
u_*=g_perp/G_DC                                             (PF-8)

whenever I-A_S is invertible and G_DC is nonzero. The shipping Kalman/S
structure contains no identity forcing G_DC=0. Thus the mean-level Fredholm
matrix is generically capable of supporting the required nonzero AW component.

The same conclusion holds for the full multi-input/multi-epoch period:
measurement corrections break the open LIN-chain integral conservation.
There is no source-uniform left annihilator of the literal accelerometer/S/
magnetic correction columns whose pairing with q_B is forced nonzero by the
gravity/field offset. Therefore PF-5 cannot be declared inconsistent from
the current state equations and boundedness assumptions alone.

### Add periodic physical moments using the SAME waveform

Periodic physical velocity and position impose exact linear moment equations
on p:

sum_k dt_k a_phys,k=0,                                      (PF-9)

and the discrete first acceleration moment equals the required velocity/
position boundary term (zero for a fully periodic translational state after
choosing the period boundary consistently). Bounded-potential adds the next
integrated moment constraint.

On PF-1, the exact accelerometer innovation is

r_a,k = a_phys,k - a_hat_w,k - b_hat_a,k + gravity/frame/lever terms. (PF-10)

Thus a zero-mean physical acceleration need not have zero-mean innovation:
the gravity-sized biased nominal AW in PF-1 contributes a fixed mean offset.
The literal accelerometer correction can therefore supply the nonzero mean AW
increment required to balance OU decay.

After eliminating the two transverse physical acceleration components with
the exact manifold-invariance equation, write

a_perp,k=f_k+g_k u_k,                                      (PF-11)

where u_k is the remaining longitudinal physical acceleration component and
f_k,g_k are generated by the SAME rho orbit. PF-9 and the higher periodic
moments become finite linear functionals of the scalar sequence u_k. Their
uncontrollable pointwise component is

chi_k=det(g_k,f_k).                                         (PF-12)

Current MARINE MOTION, gravity span and magnetic service impose no sign or
nonzero-moment law on chi_k. Consequently the periodic moment equations do not
supply a source-uniform Fredholm contradiction either. They may constrain a
particular candidate waveform, but the theorem assumptions do not make PF-5
inconsistent for every admissible waveform.

### Reinsert the delayed tuner/covariance fixed point

The remaining condition is the genuinely coupled one. Let

rho=Rho_ship[p,x]                                           (PF-13)

be the literal private-front-end, period/variance estimator, smoothed
tau/sigma_aw, SpectralMSE R_S, tau-scaled T_S, delayed commit, Riccati/Joseph,
S/sync/magnetic/event orbit. After eliminating x by PF-5, a periodic LaSalle
obstruction is exactly a solution of

rho_L=rho_0,                                                (PF-14)
q_PF in Range(M_PF(p,rho)),                                 (PF-15)
rho=Rho_ship[p,X_PF(p,rho)],                                (PF-16)
physical periodic moment/service constraints on p.          (PF-17)

This is smaller than the previous generic zero-dynamics problem. It contains
no independent residual, covariance, gain or tuner variables.

### Verdict

The focused periodic/Fredholm calculation does NOT analytically exclude the
surviving field-aligned base execution under the current assumptions.

It also does NOT construct one. What it proves is:

1. the linear MEKF mean subsystem and periodic physical moment equations do
   not contain a structural contradiction with PF-1;
2. the only unresolved self-consistency is PF-14--PF-17, especially whether
   the SAME physical waveform makes the delayed tuner/covariance orbit produce
   the gains required by PF-15;
3. therefore the nonzero pure-field-axis LaSalle candidate remains
   UNPROVED-REACHABLE, not an admissible established counterexample.

For the dissipativity architecture this is a sharp decision point. A theorem
excluding PF-14--PF-17 would make Inv{D=0} trivial and yield finite-window
strictness by DL compactness. A strict-margin solution of PF-14--PF-17 would
be a genuine shipping-reachable zero-dissipation mode and would refute strict
contraction under the present assumptions.

The next analytical step, if the LaSalle route is continued, must attack the
tuner/covariance fixed point itself rather than derive another mean-state or
kinematic moment identity. In particular, substitute the period-scaled
front-end variance and SpectralMSE law into the periodic Riccati/S scheduler
map and test whether PF-15 can be invariant under that SAME map.


## SpectralMSE tuner + Riccati/S periodic fixed-point test

This calculation substitutes the deployed tuner law into the surviving
field-aligned periodic obstruction. It corrects one misleading phrase in the
previous PF conclusion: for a prescribed physical waveform the tuner and
covariance are not a mutually coupled algebraic fixed point. Shipping is
triangular in this part of the chronology.

### 1. The measurement-only tuner is upstream of covariance

Let p be an L-periodic accepted physical IMU waveform, including the attitude
and physical acceleration seen by the private level/front-end path. Let eta
collect the private frequency smoother, period-scaled acceleration band,
variance estimator, stillness state, tau/sigma/RS EMA states, pending commit
and scheduler progress.

The literal tuner recursion is

eta_(k+1)=T_k(eta_k,p_k),                                   (TR-1)

and contains NO P or Kalman gain input. Its operating-point targets are

f_k = clamp(f_tuner,k),
tau_t,k = clamp(c_tau/(2 f_k)),                             (TR-2)
sigma_t,k = clamp(c_sigma sqrt(max(var_band,k-var_noise,0))), (TR-3)

with the documented startup/floor branches. The deployed SpectralMSE target is

RS_t,k = clamp[
 C_MSE (2 r_a)^(1/14)
 sigma_aB,k^(6/7) tau_t,k^(24/7) / sqrt(TS(tau_t,k))
],                                                         (TR-4)

where sigma_aB=sigma_t/c_sigma and

TS(tau)=clamp(c_T tau, TS_min, TS_max).                    (TR-5)

The applied EMA states obey

tauA_(k+1)   =(1-alpha_k) tauA_k+alpha_k tau_t,k,
sigmaA_(k+1) =(1-alpha_k) sigmaA_k+alpha_k sigma_t,k,
RSA_(k+1)    =(1-beta_k) RSA_k+beta_k RS_t,k,              (TR-6)

with 0<alpha_k,beta_k<=1 on every valid Live sample. A candidate is staged
after processing sample k and, when activation cadence permits, committed at
the beginning of a later sample. The delay is therefore a finite phase state
inside eta, not an algebraic same-sample loop.

For any fixed periodic target sequence on one fixed clamp/event stratum, the
scalar affine EMA monodromy has multiplier

q_tau=prod_k(1-alpha_k), q_sigma=q_tau,
q_RS=prod_k(1-beta_k),                                    (TR-7)

strictly in [0,1). Hence each applied tuner channel has a UNIQUE periodic
orbit. Explicitly, for x_(k+1)=a_k x_k+(1-a_k)t_k,

x_0^* =
 [sum_(i=0)^(L-1) (1-a_i) t_i prod_(j=i+1)^(L-1) a_j]
 /[1-prod_(j=0)^(L-1)a_j].                                 (TR-8)

The delayed commit/pending state merely shifts/samples this unique periodic
candidate sequence once the scheduler phase itself is periodic. Clamps are
nonexpansive and do not create a contradiction; on a fixed active clamp face
they replace the target by the corresponding constant boundary value.

Therefore the actual coupled relation (tau,sigma,RS,TS) is restrictive--TR-4
and TR-5 prohibit independent extrema--but it does NOT generically obstruct a
periodic tuner orbit. For a prescribed periodic physical waveform it
constructs one.

### 2. SpectralMSE makes the tuple lower-dimensional, not inconsistent

Away from cadence clamps, TS=c_T tau and TR-4 reduces exactly to

RS_t =
 C_* sigma_aB^(6/7) tau^(41/14),                           (TR-9)

because 24/7-1/2=41/14. Thus the instantaneous target tuple lies on a
two-dimensional graph parameterized by (f,var_band), and after EMA/delay the
applied tuple lies on the causal filtered image of that graph.

This is the correct replacement for a four-dimensional box. In particular a
proof may not choose a large tau with an independently small RS or unrelated
sigma. But TR-9 has no sign/equality relation involving the field-alignment
mean condition P_B a_w=P_B g. The tuner sees the physical/private-front-end
band, not the nominal AW DC offset itself.

A periodic physical waveform can have zero translational-acceleration mean
while having positive band variance and a finite frequency. TR-2--TR-9 then
produce finite positive periodic tau,sigma,RS,TS. Nothing in SpectralMSE
forces the nominal AW transverse mean to zero.

### 3. Covariance is downstream: periodic Riccati/S map

Given the periodic applied tuner/scheduler orbit eta^*(p), the covariance map
is deterministic:

P_(k+1)=R_k(P_k; p_k,eta_k^*),                              (TR-10)

where R_k is the literal sequence of prediction with Q(tau,sigma), pending AW
PSD floor, due S Joseph update using RS and TS, accepted accelerometer and
magnetic Joseph updates, reset congruence/projection, and covariance sync
chronology.

The tuner does not read P, so TR-10 cannot invalidate TR-1--TR-9 by feedback.
The periodic covariance condition is simply

P_0=R_[0:L)(P_0;p,eta^*(p)).                               (TR-11)

Every constituent covariance operation maps PSD matrices to PSD matrices.
On the retained compact covariance tube assumed by the regional proof,
R_[0:L) is continuous on each fixed event/gate stratum. Therefore Brouwer
would give a periodic covariance fixed point IF one has a convex compact
forward-invariant covariance set for this same-history map. The current proof
has not established that invariant set independently of the old covariance
work, so TR-11 existence is not promoted here.

More strongly, standard periodic Riccati intuition suggests uniqueness/
attraction when the periodic pair is stabilizable/detectable, but importing
that theorem without verifying the literal reset, pseudo-update, sync and
rank conditions would be unjustified. The correct status is:

- the tuner orbit exists uniquely for a prescribed periodic target/event
  sequence;
- the covariance periodic orbit is a downstream Riccati fixed-point problem;
- no algebraic contradiction between SpectralMSE and Riccati/S exists.

### 4. Reinsert the field-aligned mean Fredholm equation

The covariance orbit matters because it determines the gains in

M_PF(p,P) x_0=q_PF(p,P).                                   (TR-12)

Thus the COMPLETE periodic obstruction is now

eta=eta^*(p) from TR-1--TR-9,                              (TR-13)
P_0=R_[0:L)(P_0;p,eta^*(p)),                               (TR-14)
q_PF(p,P) in Range M_PF(p,P),                              (TR-15)
physical periodic moment + service constraints on p.        (TR-16)

This is substantially smaller than the previous formulation. There is no
independent tuner fixed-point unknown at all after p is chosen. The only
nonlinear internal fixed point is P (plus scheduler phase if the candidate
period does not already close it), and P enters the obstruction only through
the literal gains/coefficient matrices in TR-15.

### 5. Can TR-13--TR-16 be excluded analytically now?

No. Substitution of the actual SpectralMSE law does not create a contradiction.

Indeed it removes degrees of freedom in the favorable direction for rigor:
tau,sigma,RS,TS are a deterministic causal functional of p. But the resulting
positive finite schedule is fully compatible with a periodic Riccati/S
recursion. There is no equation in TR-4/TR-9 that conflicts with the required
gravity-sized nominal AW mean, because the tuner is measurement-only and its
variance channel removes the noise floor and tracks the physical wave band.

Conversely this is NOT a construction of a shipping counterexample. To prove
one, one must exhibit p and a strict-margin periodic P satisfying TR-14 such
that TR-15 holds and all physical/service gates remain admissible.

### 6. Important consequence for the dissipativity strategy

The hoped-for second-stage exclusion

field-aligned mean
 + SpectralMSE tuner
 + Riccati/S
 => contradiction

does not follow structurally. The tuner is triangular and contractive, so it
is not the mechanism that removes the pure field-axis LaSalle mode.

The remaining decisive object is the COMPOSED map

G_L(p,P):=
 [ R_[0:L)(P;p,eta^*(p))-P ;
   Pi_left(p,P) q_PF(p,P) ],                               (TR-17)

where Pi_left projects q_PF onto the left-null complement of M_PF (or use the
equivalent stable-monodromy residual). A genuine periodic obstruction is a
zero of G_L together with physical/service constraints.

This is now a finite-dimensional same-history fixed-point/root problem for
(p,P), with the tuner analytically eliminated. An analytical exclusion would
need a sign/degree/range theorem for G_L; a computer-assisted proof could
interval-enclose G_L over the admissible periodic waveform/covariance
parameterization. A strict-margin zero would be a genuine shipping-reachable
LaSalle obstruction.

### Verdict

The delayed SpectralMSE tuner + Riccati/S substitution does NOT exclude the
field-aligned periodic fixed point from the current assumptions. It also
shows why: the tuner is upstream and uniquely determined by the physical
waveform, while covariance is downstream. Their coupling removes independent
extrema but supplies no contradictory equality.

Therefore the dissipativity/LaSalle route has reached its genuine theorem
boundary: strictness is equivalent to excluding zeros of TR-17 over admissible
same-history physical/covariance periodic or recurrent histories. Further
scalar tuner or covariance bounds would again discard the linked structure.


## Implicit periodic Riccati sensitivity and covariance elimination

This calculation differentiates the literal periodic covariance map on one
fixed accepted-event/clamp/scheduler stratum. It does not replace shipping by
a generic DARE.

Let the tuner already be eliminated as eta=eta^*(p) by TR. Define the
one-period covariance residual

F_P(P,p):=R_L(P;p,eta^*(p))-P.                              (RS-1)

A periodic covariance orbit satisfies F_P(P^*(p),p)=0.

### Exact Frechet derivatives of literal covariance operations

For a prediction with fixed coefficient F and process covariance Q,

Psi_pred(P)=F P F'+Q,                                      (RS-2)

so

D_P Psi_pred[Delta]=F Delta F'.                            (RS-3)

Parameter/physical variations contribute
dF P F'+F P dF'+dQ in the forcing derivative D_p Psi.

For any accepted linearized measurement/Joseph update with fixed H,R,

Psi_m(P)=P-P H'(H P H'+R)^-1 H P.                          (RS-4)

Let K=P H'S^-1, S=H P H'+R, and L=I-KH. Direct differentiation, including
the derivative of S^-1, gives the exact identity

D_P Psi_m[Delta]=L Delta L'.                               (RS-5)

Thus the covariance tangent does NOT require differentiating K separately.
Variations of H,R caused by p/eta enter only the affine tangent forcing
D_p Psi_m.

A fixed reset/frame covariance congruence

Psi_G(P)=G P G'                                             (RS-6)

has D_P Psi_G[Delta]=G Delta G'. The derivative with respect to p includes
dG P G'+G P dG'.

The default pending AW covariance floor is applied inside prediction as a PSD
increment that depends on the pre-floor AW marginal and the tuner target. On a
fixed active floor branch it is an affine map in P; therefore its P derivative
is an explicit linear projection modification. The legacy immediate AW block
replacement likewise has a projection derivative (zero on the replaced AW
marginal, identity on untouched blocks, with the literal cross-block policy).
No branch may be differentiated across its switching boundary; sensitivity is
stratum-local and one-sided at a boundary.

Symmetrization is the linear projection Delta -> (Delta+Delta')/2. Numerical
PSD repair/projection is inactive on a strict-margin analytical stratum; if it
activates, smooth IFT is not applicable there and the branch must be treated
separately.

### Periodic tangent operator

Compose the exact operation derivatives in shipping order over the period.
This gives a linear operator on symmetric covariance perturbations

L_P := D_P R_L(P^*;p,eta^*(p)).                            (RS-7)

Without an active marginal-replacement/floor projection, every operation is a
congruence, so

L_P[Delta]=A_c Delta A_c',                                 (RS-8)

where A_c is the chronological product of prediction F, accepted-update
closed-loop factors (I-KH), and reset/frame G factors over the SAME period.
With active AW floor/sync projection, insert its literal linear projection
between these congruences; RS-7 remains exact.

Vectorizing the pure-congruence case,

vec L_P=(A_c kron A_c) vec Delta.                           (RS-9)

Hence

rho(L_P)=rho(A_c)^2.                                        (RS-10)

Therefore a sufficient and, in the pure-congruence case, exact criterion for
local covariance fixed-point invertibility is

1 notin spectrum(L_P).                                     (RS-11)

The stronger contraction condition rho(A_c)<1 gives
rho(L_P)<1 and the Neumann inverse

(I-L_P)^-1=sum_(n>=0) L_P^n.                               (RS-12)

This is the precise covariance-tangent condition required by the implicit
function theorem. It is NOT the old state-error O1/O2 contraction claim:
A_c is the covariance sensitivity product for one fixed periodic base
execution.

With projection branches, use the spectrum of the exact composed symmetric-
matrix operator L_P; a norm bound <1 is sufficient but not necessary.

### Eliminate P locally

Differentiate F_P(P^*(p),p)=0:

(I-L_P)[dP^*]=B_p[dp],                                     (RS-13)

where

B_p:=D_p R_L(P^*;p,eta^*(p))                               (RS-14)

is the TOTAL waveform forcing derivative. B_p includes:
- direct physical dependence of prediction/measurement/reset Jacobians;
- derivative of the private tuner orbit eta^*(p);
- derivatives of tau,sigma_aw, SpectralMSE R_S and T_S including EMA and
  delayed commit;
- derivative of Q(tau,sigma), S noise/cadence and Racc vibration inflation;
- event-time derivative only inside a fixed scheduler stratum (event changes
  are nonsmooth boundaries).

If I-L_P is invertible,

dP^*=(I-L_P)^-1 B_p[dp].                                   (RS-15)

Thus P is locally a unique C1 function P^*(p) on every strict-margin periodic
stratum satisfying RS-11.

### Physical-waveform-only obstruction

Use the stable-monodromy form of the field-aligned mean residual when
available,

h_PF(p,P)=H_B(p,P)(I-A_L(p,P))^-1 b_L(p,P)-q_B(p).          (RS-16)

More generally use any smooth left-null chart of the Fredholm residual on a
constant-rank stratum. Define

H_L(p):=h_PF(p,P^*(p)).                                     (RS-17)

Then

D H_L[p][dp]
 =D_p h_PF[dp]
  +D_P h_PF[(I-L_P)^-1 B_p[dp]].                           (RS-18)

This is the requested exact same-history sensitivity. The covariance response
is linked to the waveform through the periodic Riccati equation rather than
independently bounded.

A periodic LaSalle obstruction on this stratum must satisfy

H_L(p)=0                                                    (RS-19)

plus the physical periodic moment, MARINE MOTION, IMU BIAS, MAGNETIC SERVICE
and gate constraints. Tuner and covariance are no longer independent
variables.

### Does covariance elimination itself exclude the obstruction?

No. RS-15 is an elimination theorem, not a sign theorem. If RS-11 holds, it
makes the obstruction SMALLER and smoother but supplies no reason for H_L(p)
to be nonzero. If RS-11 fails, that likewise does not prove a pathology:
covariance may have a nonsmooth/nonunique periodic branch or the unit tangent
may be removed by an active projection/event change.

What the calculation does establish is the exact fork:

1. **Regular covariance branch:** prove RS-11 (preferably rho(L_P)<1), eliminate
   P by RS-15, and analyze H_L solely over admissible physical waveforms.
2. **Singular covariance branch:** characterize Ker(I-L_P). A periodic
   covariance perturbation in that kernel is a genuine neutral Riccati tangent;
   it must satisfy the Fredholm transversality condition
   leftKer(I-L_P) paired with B_p[dp]=0 for continuation. Treat this as a
   separate closed stratum, not by inflating independent covariance boxes.

### Important relation to the dissipativity experiment

The pure-congruence A_c in RS-8 is the same closed-loop linear factor that
appears in covariance sensitivity, but RS-10 does NOT by itself prove
homogeneous estimator contraction. Covariance tangent contraction is a local
property of the periodic Riccati orbit. Nevertheless it is exactly what is
needed to remove P from the surviving LaSalle fixed-point equations.

Therefore the next quantitative test is now concrete and much smaller:
evaluate/prove the spectral condition for L_P on a candidate strict-margin
periodic field-aligned stratum. If rho(L_P)<1 uniformly there, all covariance
degrees of freedom disappear and the theorem obstruction becomes H_L(p)=0
in physical waveform space alone. If a unit covariance tangent is forced by
the field-aligned geometry, that identifies a new genuine equality mechanism.


## Spectral test of the periodic covariance tangent on the field-aligned branch

This calculation tests RS-11 structurally on the sole surviving LaSalle
candidate. The result is negative for strict covariance-tangent contraction:
field alignment itself supplies a unit closed-loop state mode on the
pure-congruence branch.

Let the base periodic shipping execution satisfy, at every required
accelerometer epoch,

P_B a_hat_w = P_B g,                                       (CS-1)

so the nominal specific force fhat_cog is parallel to the transported body
magnetic axis b. Let

r_k=(b_k,0_bg,0_v,0_p,0_S,0_aw,0_ba)                      (CS-2)

denote the normalized pure field-axis attitude homogeneous vector in the
current local error coordinates (with the literal frame scaling understood).

### Operation-by-operation transport

1. **Prediction.** On the zero fresh-source homogeneous dynamics, the AG
   prediction/reset transport maps the field-axis attitude class into the
   next transported field-axis class:

   F_k r_k = r_k^-                                         (CS-3)

   up to the coordinate normalization/transport already used in the magnetic
   compatibility lemmas. No LIN/AW/BA component is generated on the exact
   q=0 equality branch.

2. **S pseudo-update.** H_S r=0 because r has no S/LIN component. Hence

   (I-K_S H_S) r = r.                                      (CS-4)

3. **Accelerometer update.** The homogeneous row is
   H_a r=-[fhat_cog]x b=0 by CS-1. Therefore, independently of K_a,

   (I-K_a H_a) r = r.                                      (CS-5)

   This is exact: covariance cross terms cannot create contraction when the
   measurement row itself annihilates the vector.

4. **Magnetic update.** H_m r=-[B]x b=0, so

   (I-K_m H_m) r = r.                                      (CS-6)

5. **AW covariance sync/floor.** On the pure-congruence covariance branch this
   is absent/inactive by hypothesis; mean-state AW sync does not act on the
   attitude vector. Active covariance projection branches are treated below.

6. **Quaternion/reset/frame change.** The literal reset/frame Jacobian G
   changes coordinates but transports the same physical infinitesimal
   field-axis rotation:

   G_k r_k = r_(k+).                                       (CS-7)

Thus every accepted closed-loop factor in A_c either fixes r or transports it
to the next representation of the same physical field-axis rotation.

For one periodic base orbit, the physical magnetic/reference/frame state
returns after L, so in the same root coordinates

A_c r_0 = r_0.                                              (CS-8)

Hence

1 in spectrum(A_c),  rho(A_c)>=1.                           (CS-9)

No covariance values or tuner parameters enter this conclusion beyond their
role in maintaining the field-aligned base branch and accepted event
chronology.

### Consequence for the pure-congruence covariance tangent

RS gives

L_P[Delta]=A_c Delta A_c'.                                 (CS-10)

Set Delta_0=r_0 r_0'. Then by CS-8,

L_P[Delta_0]=Delta_0.                                      (CS-11)

Therefore

1 in spectrum(L_P),                                        (CS-12)

and I-L_P is singular. The regular implicit elimination RS-15 cannot be used
on the exact field-aligned periodic equality branch.

This is a GENUINE neutral Riccati tangent associated with the same geometric
unobservability as the nonzero LaSalle mode, not an artifact of independently
boxed covariance.

### Does this imply a family of periodic covariance fixed points?

No. A unit derivative does not by itself prove a nonlinear continuum of fixed
points. Continuation of a periodic covariance solution under waveform
variation must satisfy the singular Fredholm condition

<Lambda, B_p[dp]>=0                                        (CS-13)

for every left unit tangent Lambda in Ker(I-L_P^*), and higher-order terms may
remove or bifurcate the branch. The periodic covariance itself can remain
unique even with a unit derivative at a nonhyperbolic fixed point.

Nor does CS-12 construct the base field-aligned physical execution. It is
conditional on CS-1 being shipping reachable.

### Active AW floor/sync projection

The identity CS-11 is a statement about the pure-congruence covariance
sensitivity branch. The default pending AW floor acts only on the AW covariance
sector. Since Delta_0 has support purely in the attitude field-axis sector, an
AW-only projection leaves Delta_0 unchanged. Therefore the unit tangent
survives the default AW floor/sync derivative as well, provided that operation
does not explicitly zero attitude/AW cross terms involving Delta_0 (there are
none for Delta_0=r r').

Likewise S covariance operations have H_S r=0 and their Joseph tangent fixes
Delta_0. Accelerometer and magnetic Joseph tangents fix it by CS-5--CS-6.
Thus the actual deployed covariance-floor chronology does not remove this
particular unit tangent.

A full covariance reset that explicitly overwrote the attitude marginal could
remove it, but no such recurring shipping operation exists in Live A21.

### Stronger conclusion

The proposed spectral route cannot globally eliminate covariance on the
field-aligned equality stratum:

rho(L_P)<1 is FALSE there.                                  (CS-14)

Indeed the exact neutral tangent is

Delta_0 = r_field r_field'.                                (CS-15)

This is valuable because it aligns the covariance and dissipativity pictures:
the same physical field-axis gauge is simultaneously
- a zero-dissipation homogeneous error direction;
- a unit closed-loop state multiplier;
- a unit periodic covariance tangent.

There is therefore no hidden covariance contraction capable of eliminating
the surviving LaSalle mode once exact field alignment holds.

### What remains the correct theorem question

The entire strict-stability issue is now upstream of this neutral geometry:

Can the literal shipping BASE execution satisfy CS-1 indefinitely under the
physical/front-end/tuner/mean equations?

If NO, the field-aligned stratum is unreachable and the unit tangent never
belongs to an admissible recurring execution; DL compactness can yield strict
finite-window dissipation.

If YES, strict contraction of full attitude error is impossible under the
current assumptions because the actual measurements are geometrically
collinear on that execution. One must then accept the field-axis gauge, add an
assumption excluding persistent gravity/magnetic collinearity of the NOMINAL
specific force, or weaken the theorem.

Further covariance spectral estimates cannot decide this. The next proof work
should return to reachability/exclusion of CS-1 in physical waveform space,
using the shipping-closed equation, rather than seek rho(A_c)<1 on a stratum
where an exact unit mode is now proved.


## Mahony/adaptation-proxy excitation test on the field-aligned branch

This calculation tests the proposed bridge

persistent field alignment => positive private-proxy band energy.            (PX-1)

The literal shipping code does NOT support PX-1 from the current assumptions.

The default tuner input is produced by VerticalAccelComplementary, a private
measurement-only Mahony observer. For conditioned body specific force f_B and
its private body-to-NED quaternion R_M, shipping reports

a_proxy = -((R_M f_B)_z + g).                              (PX-2)

This proxy reads no MEKF state. It is then passed through the period-scaled
AdaptiveWaveBandPass before the variance estimator used for sigma_aw.

Consider an ideal admissible rigid-body motion with zero translational
acceleration, bounded arbitrary attitude motion, exact/calibrated gyro, and
the private Mahony observer initialized at the true tilt. Then

f_B = R_true'(-g e_z),                                     (PX-3)

and the gyro propagation preserves R_M=R_true (up to irrelevant yaw gauge).
The Mahony accelerometer correction is zero because the measured gravity
direction agrees with its predicted vertical. Hence

(R_M f_B)_z=-g                                             (PX-4)

at every sample and therefore

a_proxy=0.                                                 (PX-5)

This remains true under bounded roll/pitch excitation as well as yaw:
attitude span by itself does not create gravity leakage in a correctly
gyro-propagated levelled measurement. Zero translational acceleration also
satisfies the bounded velocity/displacement/potential and translational jerk
parts of MARINE MOTION; an allowed moving attitude episode can supply the
attitude-span requirement. Thus the current physical assumptions contain no
positive lower bound on proxy energy merely from attitude excitation.

The period-scaled band is linear on a fixed frequency schedule, so PX-5 gives
zero band signal after transients. The variance channel then sees only its
explicit measurement/noise/startup floor; it does not acquire a positive
field-alignment-dependent energy. Consequently no theorem of the form

E_proxy(W)>=E_min(g_perp,Delta_R,...)>0                    (PX-6)

follows from MARINE MOTION + the internal MEKF condition
P_B a_hat_w=P_B g alone.

This does NOT construct the complete field-aligned shipping pathology. The
MEKF condition is internal and the same raw physical history must also solve
the shipping self-consistency/Fredholm equations that sustain its biased
nominal a_hat_w. The point of PX-3--PX-5 is narrower and decisive: the private
Mahony/adaptation pipeline supplies no independent implication from attitude
motion or field alignment to positive proxy band energy. Any positive proxy
lower bound must first prove that the FULL shipping self-consistency equation
forces either

(a) nonzero translational acceleration in the proxy vertical channel, or
(b) a nonzero private-Mahony tilt tracking error/leakage,

and must quantify that forcing from the same physical history. Neither is
currently proved.

Therefore it is invalid to insert a positive E_min into the tuner law before
solving that upstream reachability relation. The only unconditional tuner
lower scale on PX-5 is the implemented variance/startup floor, which is a
design floor and is not evidence against the pathology.

### Consequence for the proposed amplitude route

The desired chain

field alignment -> Mahony leakage -> sigma_aB floor -> tuner response

breaks at its first arrow under the present theorem assumptions. The correct
same-history chain remains

field alignment + literal MEKF mean recursion
 -> required raw physical history
 -> private Mahony proxy
 -> period-scaled band/variance
 -> tuner.                                                  (PX-7)

Thus the next proof step must return to the shipping-closed mean equation and
derive what RAW physical acceleration/gyro history is forced by sustaining
P_B a_hat_w=P_B g. Only after that forced history is known can PX-2 be applied
to obtain a proxy-energy bound. Treating generic attitude span as proxy
excitation would repeat the proof-relaxation error this branch was intended
to avoid.


## Step 6: field-aligned shipping recurrence -> required physical IMU -> private Mahony proxy

This calculation removes the accelerometer innovation as an independent
variable. It is local to an accepted accelerometer epoch on one literal
same-history base execution, after prediction and any due S pseudo-update and
before the accelerometer correction.

Let aS_k be the world AW mean at that instant, R_k the literal world-to-B'
rotation, ell_k the literal lever-arm term, and bA_k the temperature-corrected
accelerometer-bias mean. Shipping predicts

fhat_B,k = R_k(aS_k-g)+ell_k+bA_k.                          (PH-1)

Let K_aw,k be the three AW rows of the ACTUAL accelerometer Kalman gain
computed from the carried covariance after the same prediction/S chronology.
For conditioned/de-heeled measured specific force f_B,k, shipping updates

aA_k = aS_k+K_aw,k(f_B,k-fhat_B,k).                        (PH-2)

Quaternion injection/reset does not change the AW mean itself. Let P_B denote
the fixed-world projector transverse to the qualified geomagnetic reference.
Exact field alignment at the post-accelerometer boundary requires

P_B(aA_k-g)=0.                                              (PH-3)

Substitution gives the literal physical-measurement equation

M_k(f_B,k-fhat_B,k)=q_k,                                   (PH-4)
M_k:=P_B K_aw,k,
q_k:=P_B(g-aS_k).                                           (PH-5)

No residual has been chosen: PH-4 is an equation for the actual conditioned
body specific-force sample.

### Required physical sample

On a regular rank-two branch, rank M_k=2. Let N_k be a unit vector spanning
Ker M_k and let M_k^dagger be any fixed right inverse on Range(P_B), for
example the metric Moore-Penrose inverse in the declared physical scaling.
Then every physical sample capable of restoring field alignment is exactly

f_B,k^req =
 fhat_B,k + M_k^dagger q_k + N_k zeta_k.                   (PH-6)

The scalar zeta_k is not an invented estimator control. It is the one physical
accelerometer component that field alignment does not determine at this
instant. It must be generated by one admissible physical translational/
rotational motion and must satisfy all temporal velocity/displacement/jerk
constraints over the SAME history.

If rank M_k<2, PH-4 is solvable only when q_k is in Range(M_k); otherwise field
alignment fails immediately and the LaSalle mode acquires positive action.
The rank-deficient solvable case has a larger physical nullspace and must be
handled by the same range/nullspace formulation rather than a pseudoinverse
shortcut.

Undoing de-heeling/calibration and the IMU lever model maps PH-6 uniquely to
the corresponding raw physical accelerometer sample once the physical
attitude/rate/temperature history is fixed.

### Substitute PH-6 into the actual private Mahony proxy

Shipping VerticalAccelComplementary is independent of the MEKF. At sample k
its private Mahony state defines the body down row d_M,k and reports

a_proxy,k = -(d_M,k' f_B,k + g).                           (PH-7)

Substituting the REQUIRED field-alignment sample PH-6 gives

a_proxy,k =
 A_forced,k - c_k zeta_k,                                  (PH-8)

where

A_forced,k :=
 -[d_M,k' (fhat_B,k+M_k^dagger q_k)+g],                    (PH-9)
c_k:=d_M,k' N_k.                                            (PH-10)

This is the exact same-history bridge that was missing from PX.

### Does field alignment force nonzero Mahony proxy action?

Not pointwise in general.

If c_k is nonzero, the remaining physical direction can make the instantaneous
proxy exactly zero by

zeta_k^0=A_forced,k/c_k.                                   (PH-11)

If c_k=0, the proxy is independent of the remaining physical direction and
field alignment forces

a_proxy,k=A_forced,k.                                      (PH-12)

Therefore a source-uniform pointwise lower bound follows only on epochs where
|c_k| is sufficiently small AND |A_forced,k| is bounded away from zero, or
after temporal physical constraints prevent the cancellation sequence PH-11.

This result explains both previous observations:
- generic attitude excitation alone need not excite the private proxy;
- the MEKF field-alignment recurrence DOES constrain the physical sample, but
  leaves one physical scalar degree of freedom on a regular rank-two branch.

### The exact temporal cancellation test

Over a complete window W, collect zeta=(zeta_k). The physical velocity,
position and higher bounded-potential moment equations are linear/affine
constraints on the physical acceleration and therefore, after PH-6, have form

C_W zeta = d_W.                                             (PH-13)

Jerk/amplitude and MARINE MOTION impose additional convex/nonlinear bounds
zeta in Z_W. The private Mahony state is itself driven by PH-6, so d_M,k and
therefore A_forced,k,c_k depend causally on earlier zeta. For a fixed Mahony
trajectory the proxy vector is

a_proxy = A_forced-C_M zeta,                               (PH-14)

with C_M diagonal entries c_k. The exact minimum proxy energy compatible with
field alignment and physical moments is

E_proxy^FA(W)=
 min_{zeta in Z_W, C_W zeta=d_W}
 || B_W[A_forced-C_M zeta] ||^2,                            (PH-15)

where B_W is the LITERAL period-scaled adaptation band operator generated
causally by that same proxy/frequency history. In the full problem B_W and
A_forced,C_M are history dependent, so PH-15 denotes the coupled constrained
functional, not a fixed quadratic program.

The desired theorem E_proxy^FA>=E_min>0 is therefore equivalent to proving
that the zero/low-band cancellation sequence PH-11 cannot simultaneously
satisfy the physical moment constraints and the private Mahony recursion.

### Characterization of the surviving low-proxy orbit

A zero-proxy candidate must satisfy at every sample

d_M,k' f_B,k^req = -g,                                     (PH-16)

together with PH-4. On a regular branch this is the 3x3 physical system

[ M_k ; d_M,k' ] f_B,k =
[ M_k fhat_B,k+q_k ; -g ].                                 (PH-17)

If the augmented matrix is nonsingular, PH-17 uniquely determines the body
specific-force sample required to keep BOTH the MEKF field aligned and the
private proxy zero. This sample then drives the Mahony gyro/accelerometer
recursion, the physical translational dynamics, and the next tuner/covariance
state. Thus the surviving candidate is no longer a free zero-dynamics family:
it is the deterministic/implicit SAME-HISTORY orbit generated by PH-17 plus
the physical gyro history.

If the augmented matrix is singular, existence is the corresponding Fredholm
range condition; incompatibility at any epoch gives positive proxy action or
breaks field alignment.

This is the correct surviving obstruction to test.

### Consequence for the proof plan

Step 6 does NOT yet prove a positive proxy-energy floor. It reduces that claim
to a precise shipping reachability problem with one physical scalar eliminated.
The next calculation must propagate PH-17 through the private Mahony update
and physical velocity/displacement equations over a complete service/excitation
window.

There are now two clean outcomes:

1. PH-17 cannot remain physically admissible/periodic while satisfying MARINE
   MOTION, IMU BIAS and MAGNETIC SERVICE. Then PH-15 has a positive compact
   minimum E_min and Parts V--VIII of the dissipativity proof can proceed.

2. PH-17 has a strict-margin admissible recurrent solution. Then the Mahony
   proxy can remain zero/low while the MEKF stays field aligned; adaptation
   does not eliminate the LaSalle mode. That solution is the first genuine
   shipping-closed obstruction and should be constructed explicitly rather
   than relaxed away.

No tuner propagation is performed yet because the required positive proxy
energy has not been proved.


## Direct shipping-admissibility reduction of the final LaSalle candidate

The goal is now to exclude, not relax, the sole surviving nonzero
zero-dissipation mode. Step PH showed that field alignment plus zero private
Mahony proxy fixes the conditioned physical accelerometer sample by

[ M_k ; d_M,k' ] f_B,k =
[ M_k fhat_B,k+q_k ; -g ],                                 (SA-1)

on a regular branch, with M_k=P_B K_aw,k and
q_k=P_B(g-aS_k).

### Physical/Mahony constraints alone cannot exclude SA-1

Write the private Mahony down axis in world coordinates as
s_M=R_true d_M. Ignoring only the already declared calibration/lever defects,
zero proxy is exactly

d_M' f_B=-g
<=> s_M' a_phys = g(s_M,z-1) = -g(1-s_M,z).                (SA-2)

This sign-definite identity is useful but not contradictory. On a recurrent
physical orbit, integration by parts gives

g integral(1-s_M,z) dt
 = integral dot(s_M)' v dt                                 (SA-3)

(up to periodic boundary terms). Current MARINE MOTION bounds the right-hand
side but supplies no positive lower separation: s_M may equal the true down
axis.

Indeed the all-time STILL branch gives the decisive admissible physical
history

a_phys=0, omega=0, v=const (take v=0), p=const,             (SA-4)

with calibrated bounded biases. The private Mahony observer can be settled at
the true tilt, so s_M=e_z and SA-2 holds with equality; its levelled proxy and
period-scaled band are exactly zero after transients. This history satisfies
the existing bounded velocity/displacement/potential/jerk conditions and is
not removed by the updated MARINE MOTION assumption, which explicitly permits
complete stillness for arbitrary duration.

Therefore no theorem based on physical moments, attitude span, Mahony leakage
or proxy energy can exclude the final LaSalle candidate on the full declared
execution class. Any valid proof must also exclude a FILTER-INTERNAL
field-aligned equilibrium driven by the quiet physical history SA-4.

### Quiet-input shipping fixed point

On SA-4 the tuner/front end converges to its literal stillness/floor schedule;
call the resulting periodic scheduler/tuner sequence eta_Q (periodic because
S cadence and parameter-activation phase remain). Covariance then follows the
literal periodic quiet Riccati map and any recurrent candidate has P_Q on its
periodic orbit.

The mean system over one complete quiet scheduler period has an affine map

x_+ = A_Q(P_Q,eta_Q) x + b_Q(P_Q,eta_Q),                   (SA-5)

where b_Q contains the constant gravity/reference/bias measurement terms.
Field alignment imposes

H_FA x = q_FA                                               (SA-6)

at every required accelerometer boundary. Thus the quiet pathological
execution exists iff the exact finite system

(I-A_Q)x=b_Q,
H_FA x=q_FA                                                (SA-7)

is consistent together with the periodic covariance equation for P_Q.

Equivalently, in the stable mean-monodromy case,

F_Q(P_Q):=
H_FA (I-A_Q)^-1 b_Q-q_FA =0.                               (SA-8)

This is the strongest admissibility test presently available because every
physical, Mahony and tuner degree of freedom has disappeared. Only the literal
quiet covariance/gain orbit remains.

### Signed AW correction condition

Project SA-5 onto the fixed transverse gravity/magnetic direction n. A
field-aligned recurrent mean requires a_Q=n'a_w>=a0>0. Over one quiet period,

a_Q [1-Phi_OU,Q] =
sum_j n' Delta a_w,j^(S/acc/mag/reset),                    (SA-9)

with the exact chronological OU factors and corrections. The quiet physical
accelerometer sample has zero translational acceleration, so every
accelerometer innovation is generated by the estimator's own nominal
field-aligned force/bias/attitude mismatch.

For the pathology to persist, the TOTAL quiet closed-loop DC gain from that
mismatch into AW must have the sign that REPLENISHES OU decay. A direct AW
measurement term has the opposite stabilizing sign; only attitude/BA/LIN
cross-covariance and S-induced corrections can reverse the net sign.

Hence a sufficient exclusion theorem is the literal signed quiet-gain lemma

n' Delta a_w,total <= 0                                    (SA-10)

whenever a_Q>0 and the physical input is SA-4 (with defect margins). Then
SA-9 is impossible because its left side is strictly positive.

The current proof does NOT establish SA-10. Earlier exact chronology analysis
already showed why generic PSD arguments are insufficient: the post-S
accelerometer numerator contains signed cross-covariance terms and the S
AW correction can have either sign. But SA-10 is now required only on the
REACHABLE quiet periodic Riccati orbit P_Q, not on an arbitrary PSD covariance.

### Consequence

The attempt to prove the final trajectory inadmissible has reduced the theorem
to a concrete shipping-only question:

Does the actual quiet periodic Riccati/S orbit admit a field-aligned mean
fixed point SA-7?

If SA-10 (or directly F_Q(P_Q)!=0) is proved on that orbit, the quiet branch
is excluded. The MOVING branch must then be treated with the same signed
periodic fixed-point test on its waveform-generated tuner/covariance orbit.

If SA-7 has a strict-margin solution, then the final LaSalle mode is genuinely
shipping-admissible even at rest and strict full-attitude stability is false
under the current theorem statement.

No stronger conclusion is justified from the present assumptions. In
particular, complete stillness prevents using mandatory Mahony/proxy excitation
as the universal exclusion mechanism.


## Scope change: stability theorem for certified MARINE MOTION only

The main OU-III stability theorem is henceforth scoped to the MARINE MOTION
operating regime. Arbitrarily long physical stillness is NOT part of the
execution class for this theorem.

This is an explicit architecture assumption, not a claim about the current
front-end stillness flag. The intended deployed architecture places an
independent inertial-regime prefilter/state machine ahead of the OU-III motion
estimator:

physical IMU -> regime prefilter -> {CERTIFIED_STILL, MARINE_MOTION}.        (MM-1)

CERTIFIED_STILL is handled by a separate stationary branch/reset model outside
the theorem proved here. The OU-III MARINE_MOTION theorem applies only after
the prefilter has admitted the execution to the moving branch.

The prefilter must be designed so that a physically moving history capable of
violating the MARINE_MOTION theorem assumptions cannot be silently certified
as STILL. In particular the theorem does NOT identify "small Mahony proxy",
"small wave-band variance", or the existing frontEndStill flag with physical
stillness. A future implementation/proof of the regime prefilter must use
independent inertial evidence and conservative hysteresis. Its soundness is a
separate obligation.

### MARINE MOTION assumption used by this proof

For every all-time execution segment on which OU-III remains in the analyzed
moving branch:

1. the existing physical amplitude, rate, jerk, bias, bounded velocity,
   bounded displacement and bounded-potential assumptions hold;
2. MAGNETIC SERVICE and IMU BIAS hold as already stated;
3. there is a fixed excitation horizon T_E and span Delta_R_min>0 such that
   every complete interval [t,t+T_E] wholly contained in the MARINE_MOTION
   branch has physical attitude span at least Delta_R_min;
4. no arbitrarily long inertially quiescent interval belongs to this class:
   such an interval is required by architecture to be transferred to
   CERTIFIED_STILL instead.

Item 4 does not replace item 3. "Not certified still" is not itself sufficient
excitation; recurring attitude-span/service remains a first-class assumption.

Transition intervals that have neither a certified-still guarantee nor a full
T_E moving-excitation window are outside the recurring contraction statement.
They must be covered by a finite-duration bounded handoff/transition lemma in
the final hybrid theorem.

### Effect on the proof

The quiet shipping fixed-point obstruction SA-4--SA-10 is removed from the
admissible class of the MARINE_MOTION theorem by scope, not by pretending it
is unstable. It remains relevant to design and proof of the separate
CERTIFIED_STILL branch.

The dissipativity/LaSalle route now needs to exclude nonzero invariant
zero-dissipation trajectories only among histories satisfying the recurring
MARINE MOTION conditions above. The sole remaining equality candidate is the
field-axis attitude mode supported by a MOVING base execution satisfying
P_B(a_hat_w-g)=0 at every relevant accelerometer epoch.

Accordingly, future reachability work must not use the stationary
counterexample as a blocker. It must ask whether the literal field-alignment
physical-sample equation PH-4/PH-17 can persist on a history that satisfies
the recurring moving attitude-span condition and all same-history
tuner/covariance equations.

If no such MOVING execution exists, Inv{D=0} is trivial on the theorem class
and compactness yields finite-window strict dissipation. The separate hybrid
proof then combines: stationary-branch stability, bounded transitions, and
OU-III MARINE_MOTION contraction.


## Weeding the moving field-aligned trajectory: exact remaining obstruction

After scoping the theorem to certified MARINE MOTION, the stationary
counterexample is removed. This does NOT by itself make the surviving
field-aligned LaSalle trajectory impossible. The moving assumptions must be
used on the literal required-sample equation PH-4.

At every accepted accelerometer epoch, persistent field alignment requires

M_k(f_B,k-fhat_B,k)=q_k,                                   (WM-1)
M_k=P_B K_aw,k, q_k=P_B(g-aS_k).

On a regular rank-two branch,

f_B,k=fhat_B,k+M_k^dagger q_k+N_k zeta_k.                 (WM-2)

Convert this conditioned body specific force to true world translational
acceleration:

a_phys,k^W =
 R_true,k f_B,k + g + d_lever/cal,k.                       (WM-3)

Therefore

a_phys,k^W = F_k + U_k zeta_k,                             (WM-4)

with literal same-history
F_k=R_true,k(fhat_B,k+M_k^dagger q_k)+g+d_k and
U_k=R_true,k N_k.

This is the exact physical trajectory family capable of maintaining the
LaSalle field-alignment condition. All tuner/covariance dependence is inside
F_k,U_k through the carried shipping history.

### What recurring attitude span does and does not do

MARINE MOTION guarantees that on every complete T_E moving window there are
epochs i,j with physical attitude separation at least Delta_R_min. This rotates
R_true and therefore F,U. But it does NOT imply that the scalar-controlled
affine acceleration family WM-4 has nonzero energy, nonzero mean, or violates
bounded motion. A time-varying scalar zeta can exploit the rotating U_k
direction.

Indeed the already-derived moment equations are

sum dt_k (F_k+U_k zeta_k)=bounded velocity increment,      (WM-5)

and the first/second temporal moments give displacement/potential closure.
Attitude span supplies variation of U_k but no theorem currently separates
the forced moment vector generated by F from the reachable moment cone/span
generated by the scalar sequence zeta.

Thus the statement

Delta_R>=Delta_R_min => no solution of WM-1               (WM-6)

is FALSE as a structural implication. Proving it would require an additional
rank/separation property of the ACTUAL sequence U_k,F_k.

### Exact complete-window admissibility matrix

For a fixed same-history moving window W with n accepted accelerometer epochs,
stack the physical moment equations after WM-4. Let C_0,C_1,C_2 denote the
discrete integration operators for acceleration, velocity/displacement and
bounded-potential moments. Then

A_W zeta = -b_W + boundary terms,                          (WM-7)

where

A_W=[C_0 U; C_1 U; C_2 U],                                 (WM-8)
b_W=[C_0 F; C_1 F; C_2 F].                                 (WM-9)

Here U is the block-diagonal map zeta_k->U_k zeta_k. Amplitude and jerk impose
additional box/difference constraints.

Therefore the moving pathological trajectory is excluded on W iff the literal
forced moment vector lies outside the physically admissible reachable set:

-b_W notin A_W Z_W + B_W,                                  (WM-10)

where Z_W encodes acceleration/jerk limits and B_W the allowed finite boundary
increments. This is an exact same-history statement; it introduces no
fictitious control because zeta is precisely the physical component left
undetermined by the two field-alignment equations.

If WM-10 holds with a positive distance uniformly over every admitted
same-history W, compactness yields a positive action/proxy margin and removes
the LaSalle trajectory.

### Two-epoch test is insufficient

At two epochs the unknowns zeta_i,zeta_j are two scalars while only bounded
physical increments, not zero increments, are required. Nonzero attitude span
changes U_i,U_j but does not generically overdetermine the system. Magnetic
service constrains the homogeneous error direction and base magnetic
corrections; it does not directly constrain physical translation. Therefore
no valid two-epoch contradiction follows solely from attitude span plus
magnetic service.

This confirms the earlier PA/MI diagnosis in a cleaner moving-only setting:
the remaining obstruction is a LONG-HORIZON reachable-moment problem for the
literal one-dimensional physical channel U_k.

### Can current MARINE MOTION assumptions prove WM-10 uniformly?

No. They bound acceleration, jerk, velocity, displacement and potential and
require recurring attitude span, but impose no relation between translational
acceleration and attitude. Hence they permit correlated translation histories
whose scalar zeta sequence can, in principle, cancel the forced moments as
U_k rotates. The shipping equations determine F,U, but the present proof has
no invariant showing their forced moments lie outside the scalar-channel
reachable set.

This is not a proof that a full shipping counterexample exists: F,U depend on
the covariance/tuner orbit generated by the same physical history. It is a
proof that the CURRENT physical assumptions alone cannot weed it out by
attitude span, bounded motion, Mahony proxy or magnetic service.

### What would actually close the moving theorem

There are now only two rigorous routes.

**Route A -- prove a shipping-generated moment separation.**
Use the literal coupled tuner/covariance recursion to prove that every
field-aligned moving history has

dist(-b_W, A_W Z_W+B_W)>=delta_FA>0                        (WM-11)

on some uniform finite window. This is a property of the shipping-generated
F,U, not an added physical assumption. It would eliminate the pathology and
complete the LaSalle strictness lemma.

**Route B -- strengthen the MARINE MOTION/prefilter admission condition.**
Require the moving regime to certify an excitation condition that directly
excludes the scalar-channel degeneracy, e.g. a raw-IMU/translation-attitude
window condition whose literal consequence is WM-11. Such a condition must be
physically measurable/certifiable by the independent regime prefilter; simply
restating nominal MEKF observability is not acceptable.

Under the user's chosen architecture (STILL handled separately), Route B is
legitimate if the prefilter can conservatively refuse to admit ambiguous weak
motion to the contraction-certified MARINE MOTION branch. A third TRANSITION/
UNCERTAIN branch may remain bounded without a strict contraction claim.

### Current theorem status

The quiet obstruction is removed by scope. The moving field-aligned trajectory
is NOT YET weeded out by the existing MARINE MOTION assumptions. The exact
remaining test is WM-11. Further manipulation of generic zero dynamics or
covariance boxes cannot prove it.

The next productive calculation should either derive WM-11 from the literal
shipping-generated F,U over a complete moving service window, or formulate
the weakest raw-IMU admission statistic that certifies WM-11 and can be
implemented by the regime prefilter.


## Raw-IMU moving admission certificate derived from the scalar-channel obstruction

The goal is to replace the shipping-dependent separation WM-11 by a
conservative condition that an INDEPENDENT regime prefilter can evaluate from
conditioned accelerometer, gyro and its private Mahony tilt only.

WM-4 shows that every persistent field-aligned LaSalle execution has physical
world acceleration

a_phys,k = F_k+U_k zeta_k,                                 (RC-1)

with only one free scalar physical channel per sample. A certificate based on
attitude span or gyro Gram alone cannot exclude RC-1: zeta_k can vary at every
sample while U_k rotates. The measurable certificate must therefore test the
specific-force history itself against the class of one-dimensional
instantaneous explanations.

### Independent levelled raw-IMU coordinates

Let R_M,k be the private Mahony body-to-level frame, independent of the MEKF,
and let f_B,k be the conditioned/calibrated body specific force available to
the prefilter. Define the independently levelled linear-acceleration proxy

y_k := R_M,k f_B,k + g e_z.                                (RC-2)

The shipping vertical proxy is -e_z'y_k; RC-2 retains all three components.
Let omega_k be the conditioned gyro after only calibration/static bias
allowance available to the independent prefilter.

On a certified moving window W, stack y=(y_1,...,y_N). The prefilter also
knows dt_k and R_M,k, hence the physical integration operators in its own
level frame up to the declared Mahony tilt-error tube.

### Degenerate one-channel motion class

For a candidate unit body direction n, its independently levelled direction is

u_k(n)=R_M,k n.                                             (RC-3)

Define the class of scalar-channel acceleration histories

C_W(n):={ c_k+u_k(n) zeta_k },                              (RC-4)

where c belongs to the declared calibration/gravity/Mahony-error tube and
zeta satisfies the same physical acceleration/jerk and finite-window
velocity/displacement boundary allowances used by MARINE MOTION. Let Z_W be
that convex admissible scalar sequence set and B_W the finite boundary tube.

For fixed n, the best-fit weighted residual is

E_W(n):=
 min_(zeta in Z_W, boundary in B_W)
 sum_k w_k || y_k-u_k(n) zeta_k-c_k ||^2.                  (RC-5)

The nuisance c_k is minimized only over independently certified sensor/
Mahony/calibration defect bounds; it is NOT an arbitrary three-vector.

The raw moving certificate is

E_move(W):=min_(|n|=1) E_W(n).                             (RC-6)

This is the squared distance of the actual independently levelled raw-IMU
history from every physically admissible rotating one-dimensional
specific-force channel.

A practical implementation can discretize/branch-and-bound the sphere or use
the Gram relaxation below; the theorem uses RC-6 itself.

### Unconstrained Gram lower surrogate

Ignoring moment/jerk constraints and defect tubes makes the degeneracy class
larger, hence gives a conservative LOWER residual. For fixed n the optimal
pointwise zeta is u_k'u y_k (u_k is unit), so

E0_W(n)=sum_k w_k y_k' [I-u_k(n)u_k(n)'] y_k.              (RC-7)

Using u_k=R_M,k n,

E0_W(n)=E_y - n' G_W n,                                    (RC-8)
E_y=sum_k w_k ||y_k||^2,
G_W=sum_k w_k R_M,k' y_k y_k' R_M,k.                       (RC-9)

Therefore

min_n E0_W(n)=E_y-lambda_max(G_W).                         (RC-10)

This is the correct raw-IMU Gram statistic derived from RC-1. It is NOT
lambda_min of an attitude Gram. It measures how much of the levelled
acceleration energy cannot be explained by ANY single body-axis scalar
history.

Define

gamma_raw(W):=
 E_y-lambda_max(G_W).                                      (RC-11)

Then gamma_raw>0 certifies non-collinearity of the measured specific-force
history with every single body-fixed axis, after independent levelling.

Because RC-7 discarded physical moment/jerk restrictions, adding them can
only increase the best-fit residual:

E_move(W)>=gamma_raw(W)-epsilon_defect(W).                 (RC-12)

where epsilon_defect is the explicit enlargement due to Mahony tilt,
calibration, lever and bias uncertainty. A robust version follows from
||delta y_k||<=eps_k:
the distance (not squared distance) obeys

sqrt(E_move) >=
 sqrt(max(gamma_raw,0)) - sqrt(sum w_k eps_k^2).            (RC-13)

Thus a sufficient implementable admission test is

sqrt(gamma_raw(W)) >
 sqrt(sum w_k eps_k^2)+gamma_margin.                        (RC-14)

It is computed entirely from conditioned accelerometer plus private Mahony
attitude; gyro enters through R_M propagation and can be separately required
to satisfy the existing attitude-span/rate certificate.

### Relation to the shipping field-aligned obstruction

RC-14 excludes a PURE one-body-axis physical acceleration history. To imply
WM-11 for the actual affine family F_k+U_k zeta_k, one more bridge is required:
the shipping-forced term F_k must lie inside the independently certified
defect/boundary class c_k after subtracting the measured physical history.

In general F_k is NOT a small sensor defect; it contains the gravity-scale
field-alignment correction M_k^dagger q_k. Therefore RC-14 alone does not
logically imply WM-11.

This is decisive: there is no raw-IMU statistic, independent of the estimator,
that can distinguish an arbitrary physically admissible measured trajectory
from itself merely because the estimator internally represents it as
F+U zeta. If the theorem allows arbitrary bounded translation, an independent
prefilter cannot know whether that same raw history is compatible with the
internal pathological gain/covariance state.

Hence a prefilter-only certificate can close the theorem only if MARINE MOTION
is strengthened by a PHYSICAL excitation condition whose violation contains
every field-aligned F+U zeta history. RC-11 is one candidate physical condition
only after proving the missing bridge F in the admitted nuisance class; that
bridge is currently false at gravity scale.

### A certificate that is sufficient by construction

The weakest exact admission statistic that DOES imply WM-11 must include the
filter-independent physical moment model but also a declared physical
subspace class known a priori to contain every pathological forced term.
If such a class S_FA,k can be bounded from hardware/geometry alone, define
the raw distance

Gamma_FA(W)=
 dist( y,
       { s_k+R_M,k n zeta_k :
         s in S_FA(W), |n|=1, zeta in Z_W } ).             (RC-15)

Then Gamma_FA>=gamma_move>0 excludes the pathology provided one proves
R_M R_true' F in S_FA for every field-aligned shipping history. At present
the only source-uniform S_FA from existing assumptions is essentially the
full bounded-acceleration ball, making RC-15 vacuous.

Therefore the requested independent raw-IMU sufficient condition cannot be
made nonvacuous from the CURRENT physical assumptions without either:
(a) a new physically meaningful translation-attitude excitation premise, or
(b) using some estimator-derived quantity (gain/covariance/nominal force) in
the regime certificate.

### Recommended minimal augmentation

If independence from the main filter is mandatory, add a measurable physical
condition directly:

on every certified MARINE MOTION window, the independently levelled
acceleration history has a robust multi-axis residual

gamma_raw(W)>=gamma_min>epsilon_defect(W),                 (RC-16)

AND prove/assume that the pathological field-aligned family is
single-body-axis explainable in these coordinates within the defect tube.
The first half is implementable now; the second half is the missing theorem
and cannot be asserted from WM-4 because of F.

Alternatively allow the admission monitor to read ONLY the main filter's
published K_aw/aS (without feeding back into estimation). Then it can evaluate
the exact WM-10/WM-11 distance directly from raw IMU plus F,U. That is a
runtime proof monitor rather than an independent physical prefilter, and it
would reject precisely the pathological reachable-moment class with no new
physical assumption.

### Conclusion

The correct Gram derived from RC-1 is

gamma_raw= sum w||y||^2 -
 lambda_max(sum w R_M' y y' R_M).                          (RC-17)

It is a useful raw-IMU multi-axis excitation statistic and is strictly better
than an invented attitude Gram. But under the present broad MARINE MOTION
translation class, RC-17 alone cannot imply the exact shipping separation
WM-11 because the pathological family has a non-small affine forced term F.

This calculation prevents an invalid proof shortcut. To genuinely weed out
the trajectory, either certify WM-11 with a read-only estimator-aware monitor,
or strengthen MARINE MOTION with a physical translation-attitude condition
that makes the F bridge provable.


## Zero-dissipation substitution into the forced physical term F

This calculation tests the missing bridge needed to turn raw multidirectional
marine excitation into exclusion of the final field-aligned LaSalle mode.

Recall the exact required physical family on a regular rank-two branch:

f_B = fhat_B + M^dagger q + N zeta,                         (ZF-1)
M=P_B K_aw, q=P_B(g-aS), Ker M=span{N},                    (ZF-2)

and in world physical acceleration coordinates

a_phys = F + U zeta,                                       (ZF-3)
F=R_true(fhat_B+M^dagger q)+g+d,
U=R_true N.                                                 (ZF-4)

The desired shortcut would be to use complete-word zero dissipation D=0 to
show F=U beta+d_small, reducing every pathological physical history to one
scalar moving channel plus defects.

### Homogeneous equality does not zero the base innovation

D is the quadratic action of the HOMOGENEOUS estimator-error trajectory about
one literal base shipping execution. D=0 implies that the homogeneous
measurement variations vanish:

delta r_acc=0, delta r_S=0, delta r_mag=0,                 (ZF-5)

together with zero homogeneous fresh-source action and the previously derived
delta v=delta p=delta S=delta a_w=0 on the surviving field-axis mode.

ZF-5 does NOT imply that the BASE innovations vanish:

r_acc^base need not be 0,
r_S^base=-S^base need not be 0,
r_mag^base need not be 0.                                  (ZF-6)

Indeed a persistent field-aligned base mean generally REQUIRES nonzero base
accelerometer/S correction supply to replenish the gravity-scale AW component
lost under OU prediction.

Therefore none of the D=0 equality conditions sets q=0 in ZF-2. Here q is the
BASE departure of the pre-accelerometer AW mean from the field-aligned target,
not the homogeneous BA compatibility coordinate that also used q notation in
earlier sections. To avoid ambiguity define henceforth

q_FA,k := P_B(g-aS_k).                                     (ZF-7)

Then

q_FA = P_B K_aw r_acc^base                                 (ZF-8)

on an exact post-accelerometer field-aligned base execution. This quantity is
generically nonzero whenever OU/S chronology moves the pre-update base mean
off the target.

### Exact decomposition of the forced term

Choose the Moore-Penrose right inverse for clarity. Then
M^dagger q_FA lies in Range(M') and is orthogonal to N=Ker M. Hence

f_B-fhat_B =
 M^dagger q_FA + N zeta                                    (ZF-9)

is the orthogonal decomposition of the BASE accelerometer residual into:
- the unique minimum-norm component required to restore the two field-aligned
  coordinates; and
- the one-dimensional null component invisible to that restoration equation.

Thus the gravity-scale particular component M^dagger q_FA is, by construction,
TRANSVERSE to N unless q_FA=0. It cannot be rewritten as N beta.

After world rotation,

F = R_true fhat_B+g+d + R_true M^dagger q_FA.              (ZF-10)

The last term is orthogonal to U=R_true N in the Euclidean physical scaling.
Consequently the exact distance of the correction part of F from the scalar
channel is

dist(R_true M^dagger q_FA, span U)
 = ||M^dagger q_FA||.                                      (ZF-11)

If sigma_min(M)>0,

||M^dagger q_FA|| >= ||q_FA||/sigma_max(M),                (ZF-12)

and also
||M^dagger q_FA|| <= ||q_FA||/sigma_min(M).                (ZF-13)

Therefore zero dissipation does the OPPOSITE of the hoped-for collapse:
whenever the base pre-update field-alignment departure q_FA is nonzero, the
literal restoration requires a physical residual component transverse to the
remaining scalar channel.

### Can the nominal-force part cancel this transverse component?

Yes in the raw physical acceleration F. The term R_true fhat_B+g+d is not a
small defect on the base execution. Under field alignment after correction,
the pre-update fhat_B contains the OU/S-displaced nominal AW, BA and lever
terms. Its projection onto U^perp can cancel or reinforce
R_true M^dagger q_FA. No homogeneous zero-action identity fixes that sign.

This is why the raw Gram gamma_raw cannot yet be asserted positive merely
from q_FA!=0: the measurable physical acceleration is the SUM ZF-10, not the
correction residual alone.

### A sharper measurable object: innovation rather than acceleration

The quantity that DOES have an exact transverse decomposition is the literal
base accelerometer innovation

r_acc^base = f_B-fhat_B.                                   (ZF-14)

From ZF-9,

r_acc^base=M^dagger q_FA+N zeta.                           (ZF-15)

Therefore its energy transverse to the one-dimensional null channel is exactly

||Pi_Nperp r_acc^base||^2
 = ||M^dagger q_FA||^2.                                    (ZF-16)

Over a window,

sum w ||Pi_Nperp r_acc^base||^2
 = sum w ||M^dagger q_FA||^2.                              (ZF-17)

This is a strict linked identity, not a bound.

But r_acc^base and N depend on the main filter through fhat_B and K_aw.
Accordingly ZF-17 is available to a read-only estimator-aware stability
monitor, not to a fully independent raw-IMU prefilter.

### Relation to marine oscillation

Ordinary physical acceleration oscillation does not constrain q_FA away from
zero and does not prevent cancellation inside F. Multidirectional RAW physical
acceleration likewise cannot be identified with ZF-17 without controlling the
nominal-force term.

However, if the architecture permits a read-only monitor, the correct
excitation certificate is immediate: require recurring transverse BASE
innovation action

Gamma_innov(W):=
 sum_(k in W) w_k ||Pi_Nk_perp r_acc,k^base||^2
 =sum w_k ||M_k^dagger q_FA,k||^2
 >= gamma_innov>0.                                         (ZF-18)

This condition is exactly tied to the field-alignment restoration geometry.
It cannot be gamed by the free scalar zeta.

For the LaSalle exclusion one needs the complementary statement: an infinite
field-aligned zero-dissipation execution would require the same base
innovation/covariance chronology indefinitely. ZF-18 by itself does not make
the HOMOGENEOUS field-axis attitude observable, because H_acc r_field=0 under
exact alignment. Thus even positive base innovation action is not sufficient
to kill the geometric homogeneous mode. It only certifies that the base is
actively maintaining the pathological alignment.

### Verdict

The proposed bridge

D=0 + field alignment => F in span(U)+small defects         (ZF-19)

is FALSE.

The exact reason is now proved: D=0 zeros homogeneous actions, not base
innovations, and the field-aligned base execution generally contains a
nonzero transverse particular innovation M^dagger q_FA. This term is
orthogonal to the free scalar channel, not contained in it.

Therefore neither ordinary nor multidirectional raw physical oscillation can,
under the present assumptions, universally weed out the pathology by a
one-channel Gram argument.

This leaves the theorem with a genuine geometric issue: if the BASE nominal
specific force is exactly collinear with the magnetic field, the pure
field-axis attitude error is unobservable regardless of how rich the BASE
innovation or physical acceleration is elsewhere. To exclude that mode, the
MARINE MOTION theorem must directly guarantee recurring NONCOLLINEARITY of the
nominal/physical gravity-sensitive vector and magnetic field, or the estimator
architecture must supply an additional independent attitude reference.

A physically meaningful theorem assumption is therefore a recurring
gravity/specific-force excitation condition, for example a windowed lower
bound on

sum_(k in W) w_k
 || P_(b_k)^perp f_ref,k ||^2 >= gamma_col>0,               (ZF-20)

where f_ref must be defined from an independent physical proxy (not the
pathological MEKF nominal force) with a proved defect tube to the
accelerometer attitude vector used by the filter. Establishing that defect
tube is the next required bridge if the condition is to be raw-IMU
certifiable.


## Physical-to-nominal force tracking under the literal accelerometer update

This section formalizes the empirical chart observation that the OU-III
nominal specific force follows the physical accelerometer force. The result is
an exact local feedback identity plus a forced tracking recursion. It is then
tested against the field-aligned LaSalle candidate.

Let h_a(x) be the literal accelerometer measurement function and z_a the
conditioned physical body specific-force measurement. Define the base
measurement mismatch

e_k := z_a,k-h_a(x_k^-).                                    (FT-1)

On a fixed accepted-update stratum, linearize h_a at x_k^- with the literal
Jacobian

C_k=[J_att,k, J_bg,k, ..., R_wb,k (AW), I (BA), ...].       (FT-2)

The Kalman correction is

delta x_k=K_k e_k,
K_k=P_k^- C_k' S_k^-1,
S_k=C_k P_k^- C_k'+R_k.                                    (FT-3)

To first order, the post-correction mismatch against the SAME physical sample
is

e_k^+ = e_k-C_k K_k e_k
      = (I-C_k K_k)e_k.                                    (FT-4)

Using S=CPC'+R,

I-C K
 = I-CPC' S^-1
 = R S^-1.                                                  (FT-5)

Thus

e_k^+=R_k S_k^-1 e_k.                                      (FT-6)

In the R-weighted measurement metric this map is similar to

(I+R^-1/2 C P C' R^-1/2)^-1,                              (FT-7)

so all eigenvalues lie in (0,1], and they are strictly below one on every
measurement direction with positive predicted measurement covariance. If

lambda_min(R^-1/2 C P C' R^-1/2)>=mu_a>0,                  (FT-8)

then

||e_k^+||_(R^-1)
 <= alpha_a ||e_k||_(R^-1),
alpha_a:=1/(1+mu_a)<1.                                     (FT-9)

This is the exact linear tracking mechanism visible in the filter charts. The
nonlinear finite-angle/lever remainder is already bounded by the existing
actual-gain finite-error lemmas and can be appended as d_NL,k.

### Between samples: forced tracking recursion

The next pre-update mismatch is not e_k^+, because both the physical force and
the nominal prediction evolve. Write

e_(k+1)^-
 = A_e,k e_k^+ + Delta_phys,k - Delta_nom,k + d_NL,k,       (FT-10)

where Delta_phys is the actual conditioned physical specific-force change and
Delta_nom contains the literal OU prediction, attitude/gyro propagation,
scheduled S pseudo-update, BA temperature evolution and lever terms. No term
is independently chosen: all are generated by the same shipping history.

Combining FT-9 and FT-10,

||e_(k+1)^-|| <=
 ||A_e,k|| alpha_a ||e_k^-||
 + ||Delta_phys,k-Delta_nom,k||
 + ||d_NL,k||.                                              (FT-11)

On a compact MARINE branch, if the complete between-update coefficient has
q_e:=sup ||A_e,k|| alpha_a<1, then iteration gives the ISS tracking bound

||e_k^-||
 <= q_e^(k-k0)||e_k0^-||
 + sum_i q_e^(k-1-i)
   (||Delta_phys,i-Delta_nom,i||+||d_NL,i||).               (FT-12)

FT-12 is the rigorous statement that physical and nominal forces follow one
another. It is conditional on the source-uniform measurement-information
floor FT-8 and on q_e<1; neither follows merely from visual charts.

### Consequence for persistent nominal field alignment

Suppose the pathological base execution maintains

P_B h_cog,k^nom =0                                         (FT-13)

at every relevant post-update epoch, where h_cog^nom is the nominal CoG
specific force in world/body-equivalent coordinates. PN1 proved that the TRUE
physical specific force cannot remain within epsilon_f of the magnetic axis
for more than T_col.

At every recurring PN1-separated epoch,

|P_B e_k| >=
 |P_B f_phys,k|-defect_nom,k.                              (FT-14)

If FT-12 provides a uniform tracking tube

||e_k|| <= epsilon_track                                   (FT-15)

with

epsilon_track+epsilon_model < epsilon_sep,                 (FT-16)

where epsilon_sep is the recurring physical noncollinearity margin supplied
by PN1 on a chosen longer window, then FT-13 is impossible at that epoch.
Indeed nominal field alignment plus FT-15 would imply physical force lies
inside the same magnetic-axis tube, contradicting PN1.

This gives the desired physical-to-nominal bridge, but only after proving a
quantitative tracking tube smaller than the physical separation margin.

### Does the current theorem already provide FT-15 with the needed margin?

Not yet.

FT-6 proves strict contraction at each accepted accelerometer update only in
directions with a positive predicted measurement-covariance floor. The current
proof has lower floors on R and process sources, but it has NOT yet exported a
source-uniform positive lower bound mu_a for the full literal C P C' on every
MARINE accelerometer epoch. More importantly, FT-10 contains physical-vs-
nominal forcing increments; MARINE acceleration may change on the same time
scale as the estimator, so the ISS steady tracking radius in FT-12 need not be
small compared with the geometric gravity/magnetic separation without a
quantitative calculation.

Thus the charts are consistent with and explained by FT-6, but they are not a
proof of the source-uniform inequality FT-16.

### Stronger complete-window form

A per-sample small tube is unnecessary. Over a complete moving window define
the chronological tracking transition

Phi_e(k,j)=product_(l=j)^(k-1) A_e,l R_l S_l^-1.            (FT-17)

Then

e_k=Phi_e(k,k0)e_k0+
 sum_i Phi_e(k,i+1) w_i,                                   (FT-18)

where w_i is the literal physical-minus-nominal prediction forcing including
S/OU/BA/attitude terms and nonlinear defects.

The theorem needs only that every interval longer than T_col contains an epoch
at which

||P_B e_k|| < |P_B f_phys,k|-epsilon_model.                (FT-19)

A source-uniform complete-window gain bound can be much sharper than FT-12,
because it retains the actual coupled tuner/covariance chronology and marine
forcing instead of independently maximizing alpha_a and w_i.

This is now the correct quantitative proof target. It is a tracking-gain
problem, not an observability/kernel problem.

### Status

The literal Kalman algebra PROVES that accepted accelerometer updates contract
physical-vs-nominal specific-force mismatch according to FT-6. Bounded
physical velocity PROVES recurring physical force/magnetic noncollinearity.
Together they reduce exclusion of the pathological nominal field alignment to
one quantitative margin:

complete-window tracking error < recurring physical noncollinearity margin.
                                                                    (FT-20)

FT-20 is not yet discharged numerically/source-uniformly. If it holds on the
declared MARINE class, the field-aligned LaSalle execution is impossible,
Inv{D=0}={0}, and dissipativity+compactness yields finite-window contraction.

If FT-20 fails for the declared broad marine envelopes, the theorem needs
either a tighter certified-MARINE admission envelope or a weaker tracking
claim; it does NOT follow that the pathology is reachable.


## Quantitative complete-window tracking-margin test

The desired exclusion is

tracking mismatch at some recurring physical-separation epoch
 < physical force/magnetic separation at that epoch.        (TM-1)

This section derives the strongest source-faithful form currently available
and identifies the exact missing quantitative input.

### 1. A strict physical separation margin needs jerk/sampling, not velocity alone

PN1 showed that true specific force cannot remain inside an arbitrarily chosen
magnetic-axis tube longer than
2 V_max/(g_Bperp-epsilon). That is a residence-time theorem. By itself it does
not supply a fixed positive sampled separation epsilon_sep: a continuous
trajectory may repeatedly leave every smaller tube by an arbitrarily small
amount.

Let
d(t):=P_B f_phys(t)=P_B a_phys(t)-P_B g.                    (TM-2)

Assume the declared translational jerk bound
|dot a_phys|<=J_max. Then d is J_max-Lipschitz. Suppose on a window of length
T all sampled/continuous transverse separations obey |d(t)|<=D. Then

P_B[v(T)-v(0)]
 =T P_B g + integral d(t)dt,                                (TM-3)

so

2 V_max >= |P_B[v(T)-v(0)]|
 >= T g_Bperp - integral |d|dt
 >= T(g_Bperp-D).                                           (TM-4)

Hence for any T>2V_max/g_Bperp,

max_(t in W)|d(t)|
 >= epsilon_sep(T):=
 g_Bperp-2V_max/T >0.                                      (TM-5)

This already gives a continuous-time excursion height. If accelerometer
accepted epochs have maximum gap h_acc and d is J_max-Lipschitz, some accepted
epoch k satisfies

|d(t_k)| >= epsilon_sep(T)-J_max h_acc.                    (TM-6)

Therefore choose T such that

epsilon_phys(T):=
 g_Bperp-2V_max/T-J_max h_acc-epsilon_field/cal >0.        (TM-7)

Then every complete MARINE window of length T contains an accepted
accelerometer epoch with true transverse specific-force magnitude at least
epsilon_phys(T).

TM-7 is the quantitative physical side of the desired margin.

### 2. Exact same-history tracking operator

For the physical-vs-nominal measurement mismatch e=z-h(x), the accepted
accelerometer correction has first-order map

e^+=D_k e + d_NL,k,
D_k:=R_k S_k^-1.                                            (TM-8)

Between accepted accelerometer epochs, carry the literal prediction, due S,
attitude/gyro, BA-temperature and lever evolution. Linearizing the measurement
mismatch gives

e_(k+1)^- = A_k D_k e_k^- + w_k,                            (TM-9)

where w_k is the TOTAL same-history physical-minus-nominal forcing plus the
finite-angle/lever remainder. Define

L_k:=A_k D_k.                                               (TM-10)

Over a complete window,

e_j^-=
 Phi(j,0)e_0^-
 + sum_(i<j) Phi(j,i+1) w_i,                               (TM-11)

Phi(j,i)=L_(j-1)...L_i.                                     (TM-12)

No gain, tuner, covariance or forcing component is independently maximized in
TM-11.

### 3. Joint action bound for the forcing term

The existing actual-gain finite-error/source-factor identity gives the correct
way to bound the accumulated w_i. Express each prediction/correction forcing
through its literal covariance/noise factor B_i u_i. After chronological
transport to epoch j,

P_j =
 Phi_P(j,0) P_0 Phi_P(j,0)'
 + sum_i Phi_P(j,i+1) B_i B_i' Phi_P(j,i+1)'.              (TM-13)

Whitening by P_j shows the horizontally stacked transported source operator has
norm at most one. Therefore for any output row C_e,j mapping state/source
defects into specific-force mismatch,

|| C_e,j sum_i Phi_P B_i u_i ||
 <= sqrt(lambda_max(C_e,j P_j C_e,j'))
    sqrt(sum_i ||u_i||^2).                                 (TM-14)

For accelerometer output this covariance factor is bounded by the predicted
measurement covariance:

C_e,j P_j C_e,j' <= S_a,j                                  (TM-15)

in PSD order after including the literal measurement model blocks. Thus in
the innovation metric,

|| accumulated estimator-source tracking defect ||_(S_a,j^-1)
 <= sqrt(sum_i ||u_i||^2).                                 (TM-16)

This is the linked complete-window bound; it avoids summing per-sample gain
norms.

### 4. Why TM-16 does not yet bound BASE physical tracking error

The u_i in the finite-error/source-factor identity are MODEL/estimation
disturbance actions: process mismatch, measurement noise/nonlinear residual,
reset defects, etc. The actual MARINE physical acceleration waveform is not a
small disturbance around a known truth trajectory in the base mean recursion;
it is the measurement signal the estimator is tracking.

In TM-9, w_i contains the physical force increment
Delta_phys-Delta_nom. There is currently no theorem bounding the action of
this BASE forcing by the stochastic/process source budget sum||u_i||^2.
Doing so would assume the conclusion that the physical waveform follows the
OU model closely enough.

The tuner chooses tau,sigma_aw from that waveform, but sigma_aw is a variance
scale, not a deterministic pathwise bound on physical acceleration increments.
Therefore the existing Joseph/source-factor identity cannot by itself produce
a deterministic source-uniform epsilon_track for arbitrary admitted MARINE
waveforms.

### 5. Exact margin condition and status

At the guaranteed physical-separation epoch j from TM-7, persistent nominal
field alignment implies

|P_B e_j| >= epsilon_phys(T)-epsilon_nominal_model.         (TM-17)

The tracking recursion TM-11 excludes the pathology if one can prove

sup_same-history
 |P_B[Phi(j,0)e_0 + sum Phi(j,i+1)w_i]|
 <
 epsilon_phys(T)-epsilon_nominal_model                     (TM-18)

for at least one such j in every T-window.

The literal Kalman correction supplies contraction D_k=R S^-1 and the source-
factor identity tightly bounds estimator/model disturbances. What remains
unbounded in TM-18 is the deterministic BASE physical forcing component of
w_i under the current MARINE MOTION class.

Consequently the quantitative tracking margin (5) is NOT proved from the
present assumptions. The obstacle is not covariance/gain independence; it is
that MARINE MOTION currently permits arbitrary bounded/jerk-limited
translational acceleration, whereas deterministic tracking error depends on
its bandwidth/amplitude relative to the estimator closed-loop bandwidth.

### 6. Minimal condition that would close it

A physically natural certified-MARINE condition is a deterministic
tracking-band envelope on the independent raw-IMU force, for example

sum_(i in W)
 || Delta f_phys,i ||^2 / q_track,i <= E_track             (TM-19)

or an equivalent Lipschitz/bandwidth bound tied to the already measured wave
period, with q_track generated from the SAME tuner schedule. If TM-19 is
chosen so that the induced complete-window gain in TM-11 is below the right
side of TM-18, nominal field alignment is impossible.

This is not an artificial observability assumption: it states quantitatively
that the marine force waveform lies inside the bandwidth/amplitude envelope
the filter is designed and tuned to track. The filter studies empirically show
this for the tested wave spectra; a theorem needs it as a deterministic
admission envelope or must derive it from a physical wave model.

Alternatively, prove a deterministic relation between the period-scaled
front-end band statistic and the pathwise forcing action in TM-19. That would
turn the existing tuner measurement into the needed certificate without a new
external assumption.

### Conclusion

The complete-window same-history calculation is now closed as far as the
current deterministic assumptions permit:

- physical separation margin: explicit TM-7;
- literal tracking transition: exact TM-8--TM-12;
- estimator/model disturbance accumulation: linked TM-14--TM-16;
- unresolved term: deterministic physical waveform forcing in TM-11.

Thus the charts' physical/nominal tracking can be proved uniformly only after
connecting the admitted marine waveform's deterministic increment energy to
the same tuner-generated tracking bandwidth. That is the next mathematical
bridge; further covariance or compatibility analysis cannot replace it.


## Literal front-end test of the deterministic tracking-action bridge

TM identified the desired bridge

sum ||Delta f_phys||_(Q_track^-1)^2
 <= C_track E_frontend(W).                                 (FB-1)

This section substitutes the actual shipping front end:
VerticalAccelComplementary -> AdaptiveWaveBandPass ->
SeaStateAutoTuner variance/frequency.

### Exact signals seen by adaptation

VerticalAccelComplementary is measurement-only. With private Mahony
body-to-NED rotation R_M and conditioned body specific force f_B, its scalar
output is

a_V=-(e_z' R_M f_B+g).                                     (FB-2)

Thus the adaptation path discards the two independently levelled horizontal
specific-force components before any band or variance operation.

AdaptiveWaveBandPass is the scalar time-varying linear recursion

l_+=q_l l+alpha_l a_V,
b_+=q_h b+alpha_h q_l(a_V-l),                              (FB-3)

with low/high corners [0.5,4] times the lagged tuning frequency, subject to
absolute/Nyquist clamps. SeaStateAutoTuner then forms debiased exponentially
weighted first and second moments of b and reports

var_B=E_w[b^2]-E_w[b]^2.                                   (FB-4)

The variance horizon is K periods of the supplied wave frequency, clamped in
seconds. The operating-point law subtracts the exact propagated bench-noise
variance of FB-3 and uses

sigma_aw,target=c_sigma sqrt(max(var_B-var_noise,0))        (FB-5)

(up to startup floor/clamp), while tau comes from the measurement-only period
estimate. SpectralMSE R_S and T_S are downstream of this same tuple.

### FB-1 is impossible on the current MARINE class

The map from the three-component conditioned physical force history to
E_frontend has a nontrivial exact nullspace even with perfect Mahony tracking.

Take a level-frame physical translational acceleration

a_phys^M(t)=A sin(omega t) e_x,                             (FB-6)

with zero vertical component, omega in any admissible marine band and A small
enough to satisfy acceleration/jerk/velocity/displacement bounds. With exact
private tilt and no vertical translation,

a_V(t)=0.                                                   (FB-7)

After the finite filter transient,

b(t)=0, var_B=0                                             (FB-8)

apart from the explicit bench/startup floors, while

Delta f_phys !=0                                            (FB-9)

and its deterministic increment/action over a nontrivial window is positive.
Amplitude A may be varied within the physical envelope without changing the
ideal frontend wave variance.

A separate admissible roll/pitch history can supply the theorem's recurring
attitude-span condition; correct gyro propagation of the private Mahony tilt
does not turn FB-6 into vertical proxy energy. Thus the counterexample is not
removed by MARINE MOTION attitude excitation.

Consequently there is no finite source-uniform C_track for FB-1 on the
current three-dimensional MARINE MOTION class.

### Band-pass nullspace makes the scalar problem noncoercive too

Even if FB-1 were restricted to vertical force, the adaptive band statistic
cannot control arbitrary deterministic vertical increment energy without a
spectral-envelope assumption. FB-3 is a high-pass followed by a low-pass.
Constant/very-slow vertical components are attenuated by the high-pass and
sufficiently fast components by the low-pass, while deterministic increment
energy may remain nonzero. The EW central variance FB-4 additionally removes
the band-output mean.

Therefore an inequality from FULL vertical pathwise increment energy to
period-scaled band variance requires the physical waveform to be restricted
to a frequency class on which the transfer magnitude of FB-3 has a positive
lower bound. The current bounded acceleration/jerk assumptions do not provide
that spectral support condition.

### Does the field-aligned pathology lie automatically in the measured channel?

No. The restoration condition is

M r_acc=q_FA, M=P_B K_aw,                                  (FB-10)

where P_B is transverse to the geomagnetic field. Its required innovation
component M^dagger q_FA is a general three-component body force. Magnetic
transverse directions are not the Mahony vertical direction. No shipping
identity proved so far forces the field-alignment-supporting physical
mismatch to have a nonzero vertical projection, let alone to lie inside
[0.5,4] f_tune.

Hence the frontend nullspace counterexample is relevant to the final
obstruction, not merely to unrelated physical forcing.

### What the tuner DOES certify

For the scalar proxy component that actually lies in the adaptive band, the
statistics are exact and useful. Conditional on a fixed coefficient history,
FB-3 is a stable linear operator B_W. The debiased EMA variance is a positive
quadratic form in its output after removal of the weighted mean. On any
declared spectral subspace X_W for which

||B_W x||_w >= beta_B ||x||_w, beta_B>0,                   (FB-11)

the measured variance supplies a coercive bound on the zero-mean scalar
component:

||x-mean_w x||_w^2
 <= beta_B^-2 E_var(W).                                    (FB-12)

The tuner then maps that SAME measured E_var causally into sigma_aw, tau,
T_S and R_S. This is a valid linked tracking-budget theorem, but only for the
measured vertical in-band subspace.

### Consequence for the LaSalle strategy

The proposed closure

arbitrary MARINE physical forcing
 -> current vertical frontend statistic
 -> uniform 3D tracking tube
 -> nominal/physical noncollinearity bridge                 (FB-13)

is mathematically false.

There are three rigorous repair options:

1. **Strengthen the certified MARINE class.** Require the physical
   acceleration relevant to the theorem to lie in a declared wave spectral/
   directional envelope with a positive projection into the private vertical
   band. Then prove beta_B>0 for that envelope and use FB-12.

2. **Make adaptation/monitoring 3-D.** Feed independent Mahony-levelled
   three-axis acceleration into corresponding period-scaled bands/statistics.
   A 3-D band-energy statistic removes the horizontal nullspace. A spectral
   support premise is still required to control out-of-band deterministic
   forcing, but it can be matched to the marine-wave class rather than to one
   vertical component.

3. **Do not use tuner variance to prove tracking.** Retain TM-18 as a direct
   same-history condition and certify it with a read-only runtime stability
   monitor using the actual physical/nominal mismatch and gains.

For the current shipping estimator and current broad MARINE assumptions,
FB-1 cannot be used to close Inv{D=0}={0}. This is a structural limitation of
the information supplied to adaptation, not a missing covariance inequality.

### Minimal theorem-compatible path

If the intended theorem is genuinely for wave-driven marine motion rather
than arbitrary bounded 3-D translation, formalize that physical class
explicitly. A self-similar marine-wave spectral envelope already motivates
the SpectralMSE tuner. If the theorem assumes that the acceleration component
relevant to force tracking has support in
[f_low,f_high]=[0.5,4] f_tune (with leakage margin) and a known nonzero
projection into the Mahony vertical channel, then the minimum transfer gain of
FB-3 on that compact normalized band is positive. FB-11--FB-12 become
quantitative and the SAME sigma_a,B used by shipping supplies the deterministic
forcing-action bound needed by TM.

Without those physical spectral/directional premises, the desired proof would
claim more than the front end observes.


## Quantitative radius-local exclusion using the declared MARINE/MAGNETIC bounds

The previous frontend detour enlarged the theorem class unnecessarily.  For
the regional theorem use the already proved retained-ball AW error estimate

|a_hat-a| <= 4 r                                            (QR-1)

(and |e_ba|<=r/40 where the force comparison needs BA explicitly), together
with the exact MARINE MOTION and MAGNETIC SERVICE constants.

Let e_B be the unit direction of
d_0=P_(b_M)^perp g_0.  The declared inclination bound gives

g_Bperp >= g_min cos(80 deg).                               (QR-2)

Let eps_F collect ONLY the already declared fixed physical/reference defects
in the scalar e_B force comparison: local-gravity variation eps_g,
near-constant-field projector variation from eps_B, lever/calibration terms,
and the exact-compatibility tolerance.  Let eps_r(r) collect radius-local
attitude/BA terms, with eps_r(r)->0 and including the r/40 BA contribution
when applicable.

If the nominal field-aligned zero-action pathology persists at every applied
accelerometer epoch t_k, QR-1 implies

e_B' a(t_k) >=
 g_Bperp - [4r+eps_r(r)+eps_F].                             (QR-3)

(The sign convention is chosen so e_B'g=g_Bperp.)

Between applied accelerometer epochs, |a_dot|<=J_max.  For a gap D_k,

integral_(t_k)^(t_(k+1)) e_B'a(t) dt
 >= [g_Bperp-eps_col] D_k - J_max D_k^2/4,                 (QR-4)

where eps_col=4r+eps_r(r)+eps_F.  Summing gaps over a window of length L and
using |v_end-v_start|<=2 V_max gives the NECESSARY condition for persistence

eps_col >=
 g_Bperp - 2 V_max/L
 - (J_max/4) [sum D_k^2/sum D_k].                          (QR-5)

For the regular 25-Hz applied cadence used by the proof,
D_k<=h_acc=0.04 s, hence

sum D_k^2/sum D_k <= h_acc                                 (QR-6)

and a persistent pathology requires

4r+eps_r(r)+eps_F >=
 M(L):=g_Bperp-2 V_max/L-J_max h_acc/4.                    (QR-7)

Therefore it is EXCLUDED whenever

4r+eps_r(r)+eps_F < M(L).                                  (QR-8)

This is the requested explicit radius inequality.

### Literal numerical substitution

The current deterministic certification table gives

V_max=5.50 m/s,
J_max=100 m/s^3.                                           (QR-9)

The magnetic theorem gives

g_Bperp>=g_min cos80deg
       =0.1736481777 g_min.                                (QR-10)

Using g_min=9.80665 m/s^2 for the nominal numerical audit gives
g_Bperp>=1.7029069 m/s^2.  This numerical use is an audit value; the theorem
retains g_min and eps_g symbolically until their certified envelopes are
instantiated.

At h_acc=0.04 s,

J_max h_acc/4=1.000000 m/s^2.                              (QR-11)

Thus, before eps_F/eps_r deductions,

L=16 s:
 M_16=1.7029069-11/16-1
     =0.0154069 m/s^2,
 r_FA,ideal(16)=M_16/4
     =0.0038517.                                            (QR-12)

L=64 s:
 M_64=1.7029069-11/64-1
     =0.5310319 m/s^2,
 r_FA,ideal(64)=0.1327580.                                 (QR-13)

L=100 s:
 M_100=1.7029069-11/100-1
      =0.5929069 m/s^2,
 r_FA,ideal(100)=0.1482267.                                (QR-14)

As L->infinity,

M_inf=1.7029069-1=0.7029069,
r_FA,ideal(inf)=0.1757267.                                 (QR-15)

The 16-s route has essentially no defect headroom and should NOT be used for
the final certificate.  The existing 100-s superword has about
0.593 m/s^2 of zero-radius headroom before declared defects and is the natural
window for this exclusion.

### Symbolic certified radius with the actual field envelopes retained

Because eps_g and eps_B are deliberately still symbolic in the theorem
contract, a fully numerical r_FA cannot honestly be emitted yet.  The
certified 100-s condition is

4r+eps_r(r)
 <
 0.1736481777 g_min
 -0.11
 -1.00
 -eps_F,                                                    (QR-16)

or

4r+eps_r(r)
 <
 0.1736481777 g_min -1.11-eps_F.                           (QR-17)

Here eps_F must be instantiated from the declared eps_g/eps_B projector
perturbation plus lever/calibration/compatibility constants.  A positive
radius exists iff

eps_F < 0.1736481777 g_min-1.11.                           (QR-18)

At nominal standard gravity the right side is about

0.5929069 m/s^2.                                           (QR-19)

Once eps_F is supplied, define r_FA as the largest nonnegative root of

4r+eps_r(r)+eps_F
 =0.1736481777 g_min-1.11.                                 (QR-20)

If the only remaining radius-dependent force defect is the BA storage term
r/40, a conservative closed form is

r_FA >=
 [0.1736481777 g_min-1.11-eps_F]/(4+1/40),                 (QR-21)

provided the numerator is positive.  At zero fixed defects and nominal
gravity this gives

r_FA >=0.5929069/4.025
     =0.147306.                                             (QR-22)

Any additional attitude-to-force term eps_att(r) must be added to the
denominator/modulus rather than omitted.

### Relation to the older jerk-collinearity audit

The existing world-frame lemma obtained, for h=1/5 and L=16 s,

sum D_k^2/sum D_k >=0.0509532 s

if physical force is exactly collinear at every correction.  A regular
25-Hz cadence has D=0.04 s, so exact collinearity is already impossible.
QR-5 is its radius-local tube version.  The very small M_16 above explains
why the exact-collinearity statement was easy while a useful finite-radius
16-s tube has little margin.  Extending to the already contemplated 100-s
superword converts the same physical mechanism into useful radius headroom.

### Analytical consequence

Fix any certified field/reference envelopes satisfying QR-18 and choose
r<r_FA from QR-20.  An infinite zero-dissipation MARINE trajectory in the
retained ball would, by the complete-word equality characterization, maintain
nominal field alignment at every applied accelerometer epoch.  QR-1 then
imposes QR-3.  But summing the literal MARINE jerk/velocity inequalities over
each 100-s block contradicts QR-8.

Hence, conditionally on the already declared field/reference envelopes,

Inv_MARINE({D=0} intersect {V<=r^2})={0}.                  (QR-23)

Compactness of the retained same-history class then yields a finite block
length m and eta_D(r)>0 with strict accumulated dissipation, hence homogeneous
block contraction.  The practical/ISS theorem still needs the existing
nonlinear/source-supply retention step, but the alleged moving
field-aligned trajectory is analytically excluded inside this radius.

### Remaining data obligation

The calculation found the actual bottleneck: not tuner variance, but the
still-symbolic certified values of eps_g and eps_B and any nonzero
lever/calibration/compatibility force defects.  They must be instantiated or
bounded tightly enough that their total eps_F is <0.5929069 m/s^2 (at nominal
gravity) on the 100-s certificate.  This is a large margin relative to normal
field/reference perturbations, but the proof must insert the declared values
rather than assume them zero.
