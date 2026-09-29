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
- **Invariance.** At the `10^-3` tilt premise `kappa_nu` is 195.8–461.7
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

The weighted least-squares identity now gives the leverage inequality

`X_g S_g^-1 X_g' <= I`.                                    (LE-4)

Consequently every occurrence of the potentially large signed history
`C_k` that lies in the observed gyro column is charged through `S_g^-1`,
not through `sup|C_k|`.  Writing the AG-process coefficient as observed
gyro-column part plus the within-step remainder,

`mathcal C = mathcal C_obs + mathcal R_D`,

the canonical reader action obeys

`M_C <= 2 F G_mu^-1 F'
       +2 F G_mu^-1 X' mathcal R_D mathcal R_D' X
            G_mu^-1 F'`.                                   (LE-5)

The first term is controlled directly by the quotient information:
`F G_mu^-1 F'<=||F||^2/s(c,r)^2 I`.  The second term contains only the
within-step defect `R_D`, not the accumulated chronological `C_k`.
The implemented gyro invariant already gives
`||R_k^-1D_k-h_k I||<=h_k(theta_k/2+theta_k^2/3)`; hence with
`theta_k<=theta_max<.007`,

`||R_D,k||<=h_k e_D`,
`e_D=theta_max/2+theta_max^2/3`.                            (LE-6)

Using (LE-4) once more on the defect-weighted rows gives the finite bound

`M_C <= 2 ||F||^2/s(c,r)^2
        [1+e_D^2 T h_max / q_T] I`.                        (LE-7)

The scalar next-kernel reader has the identical estimate with
`||F||^2` replaced by `|f_d|^2`:

`m_d,C <= 2 |f_d|^2/s(c,r)^2
          [1+e_D^2 T h_max / q_T]`.                         (LE-8)

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

Consequently the leverage identity

`X_g S_g^-1 X_g'<=I`

and the reader action use one common whitening and anchor convention.  The
literal reset transport is retained exactly; finite-series prediction error
is isolated in `R_D`; S-chain shared-source correlations are retained in
`Sigma_z`.

**Result.**  The bookkeeping identification required by the signed
`C_k`-energy bound passes in real arithmetic.  The AG-process family is
therefore bounded conditionally on the already stated positive local-tube
Schur floor `q_T(c,r)>0` and quotient floor `s(c,r)>0`.  This does not
prove those radius-dependent floors numerically; it removes the separate
whitening/anchor obstruction.

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

## Qualitative complete-word coercivity and quantitative gap

Fix a candidate retained radius r and kernel ceiling c>0, and append kernel
precision mu=1/c. The complete-word augmented joint action is continuous and
nonnegative on the normalized radius-local word/kernel class.

The zero-action results now imply trivial augmented nullspace: zero fresh
action rigidifies nuisance trajectories; four S rows kill the homogeneous
LIN/AW root; magnetic service and chronological transport remove gyro-bias
directions except the magnetic-compatible attitude line; all accelerometer
rows plus one decaying BA root leave dimension at most one; that line belongs
to the compact kernel family N(r); and the appended rank-one precision removes
it.

Therefore compactness gives a positive minimum normalized action

`g_MW(c,r)>0`,

and hence a uniform augmented information floor

`J_MW,mu(W)>=g_MW(c,r) I`.

Together with the finite known-root terminal covariance/action bounds this
proves the existence statements

`K_MW(c,r)<infinity`, `D(c,r)<infinity`.

This is qualitative coercivity, not yet a numerical certificate. Two
compactness steps remain nonconstructive: the four-S lemma proves a positive
determinant minimum but has no explicit numerical lower singular value, and
the multi-epoch magnetic/accelerometer compatibility argument proves
dimension <=1 but gives no explicit lower principal angle for the next
independent constraint. Without those two moduli, g_MW cannot be evaluated,
so neither K_MW nor D nor the contraction margin can be evaluated.

The next quantitative task is therefore an analytic/interval enclosure of
the full normalized complete-word action supplying (a) a four-S singular
value floor and (b) a minimum principal-angle floor for the complete
magnetic+accelerometer compatibility matrix over the compact retained class.
Carried-word singular values are not substitutes for these moduli.
## Explicit four-S floor and principal-angle obstruction

Extend the regular three-S selector one row backward. Four selected S times
have consecutive spacings in [8,8.156] s. Set t0=0, so t3<=24.468 s.
For the basis {1,t,t^2,psi_tau(t)}, generalized Vandermonde gives

`|det V4|=prod_(i<j)(t_j-t_i) psi_tau'''(xi)/3!`.

Since psi'''=exp(-t/tau), tau<=12, and pair distances are at least
8,8,8,16,16,24 s,

`|det V4| >= [8^3 16^2 24/6] exp(-24.468/12) = delta_det`,

with delta_det>6.82e4 in this unscaled basis. This is an analytic theorem
constant. If M4 bounds ||V4||_2, then

`sigma_min(V4)>=delta_det/M4^3`.

A fully explicit choice follows from Frobenius norm and
`|psi_tau(t)|<=t^3/6`:

`M4^2<=4[1+24.468^2+24.468^4+(24.468^3/6)^2]`.

Thus the four-S singular-value modulus is explicit.

For the full magnetic+accelerometer compatibility matrix, let C_MA stack
the projectors transverse to all pulled-back magnetic lines and the
BA-eliminated accelerometer compatibility rows. The qualitative proof gives
pointwise nullity at most one, but a uniform principal-angle floor would need

`inf_W sigma_min^+(C_MA(W))>0`.

That does not follow from current assumptions. MAGNETIC SERVICE controls
service information/gaps, not transversality between pulled-back magnetic
lines and accelerometer compatibility rows. MARINE MOTION controls physical
attitude span, while the latter rows contain nominal specific force. Admitted
collinear/sync-locked histories can approach compatibility alignment
continuously. The compact closure can therefore approach a rank-loss limit
while every nearby word still has nullity at most one.

Pointwise rank plus compactness is insufficient because sigma_min^+ is not
continuous through rank loss. Hence no positive source-uniform principal
angle is currently proved. Carried-word minima cannot fill this theorem gap.

The four-S quantitative modulus is closed; complete-word quantitative
coercivity remains blocked only by this transversality issue. A different
kernel-bounded quantity that stays regular as the one-dimensional
compatibility line rotates is required unless a new physical transversality
assumption is introduced.
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

## Explicit coupled coefficient system

For a candidate invariant pair `(r,c)`, let `B_q(c,r)` bound covariance
on the quotient of the physical kernel.  Together with
`nu'Pnu<=c` and the proved nuisance bounds, block Cauchy--Schwarz gives a
radius-local Euclidean ceiling `P<=C(c,r)I`.  Hence
`|theta|<=sqrt(C(c,r))r`, while the proved marginals give
`|e_aw|<=4r` and `|e_ba|<=r/40`.

Parameterize the two radius-dependent geometry terms by proved reader/source
constants:
`m_perp(c,r)<=A0+A1 sqrt(C(c,r)) r`, and
`delta_Q(c,r)<=Q0+Q1 sqrt(C(c,r)) r+Q2 C(c,r)r^2+Q3 C(c,r)^(3/2)r^3`.
Insert these in the local-tube G0 formulas.  With
`mu=sigma_w-m_perp/g` and `gamma=mu_+^2/(u1^2+1)`, let `s(c,r)` be the
resulting quotient singular floor.  If `R_q(c,r)` is the transported
quotient-reader action and `R_d(c,r)` the action of the same coefficients
on the terminal kernel functional, take

`K(c,r)=R_q(c,r)/s(c,r)^2`,  `D(c,r)=R_d(c,r)`.

The finite-error composition has polynomial form

`E(c,r)<=e0
 +Racc[(Fmax/2)C(c,r)r^2+4 sqrt(C(c,r))r^2]
 +Rmag[(Bmax/2)C(c,r)r^2]
 +Rreset[z1 sqrt(C(c,r))r+z2 C(c,r)r^2+z3 C(c,r)^(3/2)r^3]`.

Projection contributes zero whenever the proved projection-inactive prefix
guard applies.  The coupled certification problem is therefore exactly

`D(c,r)<=c`,
`[1-sqrt(1-1/K(c,r))]r>E(c,r)`.

A positive solution certifies
`1-rho_0(c,r)>=1/K(c,r)>0` on that invariant region.  This is currently a
symbolic reduction, not a numerical certificate: finite source-uniform
values for `A0,A1,Q0..Q3,R_q,R_d` and the associated action constants
remain OPEN.  Carried-word fitted values must not be substituted for them.
No positive `(r,c)` or numerical `rho_0` is claimed until those constants
are proved.

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
