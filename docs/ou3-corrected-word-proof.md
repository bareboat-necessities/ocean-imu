# Corrected-word covariance reduction and finite-error composition

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
action bound. Euclidean BA projection retains its separate sector/defect;
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
