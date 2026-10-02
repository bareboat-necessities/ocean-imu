# OU-III signed temporal continuation

Scope: apply moving excitation only to complete windows contained in a physical
moving episode, as defined in `ou3-regime-design.md`. No physical-STILL detector
is inferred. The smooth sin-cubed rest/motion/rest extension preserves the
quiet-input obstruction across transitions. `ou3-moving-six-pivots.md` adds a
conditional chronological transport/actual-row budget for all six pivots;
the source-uniform budget and the two-temporal-margin implication remain OPEN.
Removing indefinite rest does not improve the failed bounds on the already
excited moving diagnostic below.

The chronological continuation now also uses C=A^-1 B: literal reset inverses
have norm at most one, and resets leave C unchanged while altering later gyro
injections. The conditional prefix budget does not certify a signed temporal
margin. An exact relaxed reset sequence gives B=diag(0,0,28/625), without
proving shipping reachability or loss of all historical rows. The exact quiet
nominal subcase now has uniform historical covariance action and qualitative
homogeneous linear loss; physical stationary robustness remains open.
Actual acc/mag row groups within one prediction cell now factor with E=0.
Their geometry depends on the field pulled back through the local resets;
the signed physical-to-nominal transfer of a positive margin is still OPEN.
`ou3-world-frame-rows.md` shows that the historical rows are attitude-free in
world coordinates. The body-frame rotation/reference action whose triangle
budget fails below therefore does not enter the aggregate attitude columns;
they need only the nominal signed AW mean. Same-cell geometry also depends
on the applied magnetic cadence.


This continues the existing construction-linked proof. Each lemma below is
subordinate to `V_(j+1)<=rho V_j+c_d||d||^2, rho<1`: span and signed-balance
lemmas target the nominal temporal margins; a six-pivot lemma then supplies a
conditional historical-reader bound for the complete matrix process comparison.
None of these lemmas by itself supplies contraction, capture or retention.

## Current controlling physical supply (SLOW + FAST)

All physical b_a/b_g in the retained bias-state algebra now mean the SLOW
components. FAST errors remain the signed calibrated measurement-row forcing.
Use `ou3-imu-two-timescale.md` SF1--SF9 for the admissible domain. Both numerical
fast profiles are OPEN; neither this note's carried diagnostics nor old
amplitude-only sums establish them. Historical O1/O2/reader proposals below
are algebraic research history, not the controlling proof architecture.

For any actual row weights A_k, slow action is exactly
`(sum A_k)b_s,0 + sum_j(sum_(k>=j)A_k)w_s,j`, with linked increments and ALL
partial amplitude constraints. Fast action is exactly
`q_last(U_end-U_start)+sum_j(q_(j-1)-q_j)(U_j-U_start)`, q_k=A_k/dt_k.
Keep the signed endpoint and varying gyro/attitude/frame gains. Their norms
may be bounded using K* only after temporal qualification. Form the linked
chi_gamma on these SAME histories; a failed norm-sum diagnostic is not a
reason to fit H,C or to assert a source-uniform inequality.

## Physical span and the literal signed relations

Write Q(t) for true world-to-body rotation, R_i for the estimator rotation,
and d(t)=Q(t)e_z. Angular distance satisfies
`dist(d(s),d(t))<=Omega_max |s-t|`. Choose points attaining the diameter of
the complete physical window and a nearest sample to each. If the sample
coverage distance is eta, the spherical triangle inequality proves

`Delta_g(samples)>=alpha:=max(0,theta_E-2 Omega_max eta)`.

Some sampled pair therefore has chord at least `2 sin(alpha/2)`. In-window
samples give eta<=h_max when both ends are covered; eta<=h_max/2 requires the
corresponding endpoint or padded-cell coverage. T_E and theta_E remain
symbolic. This estimate is not strictly positive when theta_E<=2 Omega_max eta.
It never identifies nominal attitude with true attitude.

At actual accelerometer events retain exactly

`r_i^a=f_i^m-R_i(ahat_i-g e_z)-bhat_a,i`,
`f_i^m=Q_i(a_i-g e_z)+b_a,i+n_a,i`.

Writing fhat_i=R_i(ahat_i-g e_z), subtraction on the same physical history gives

`g(d_j-d_i)=Delta(Qa)-Delta fhat+Delta b_a-Delta bhat_a+Delta n_a-Delta r^a`.

Thus the sampled pair satisfies

`||Delta fhat+Delta bhat_a+Delta r^a-Delta(Qa)||
 >=2g sin(alpha/2)-min(2B_a,D_a,s |t_j-t_i|)-2N_a`.

This constrains a signed combination. It does not lower bound nominal
force/field angle: physical acceleration, actual bias updates and innovations
cannot be deleted. The magnetic relation remains

`r_i^m=m_i^m-R_i B_i`, `m_i^m=Q_i B_i^true+b_m,i+n_m,i`.

B_i is the actual nominal reference, including adaptation/refinement and hard
events. MAGNETIC SERVICE remains the actual normalized heading/axial-bias
information sum with its preceding corrected transport, sensitivity and
innovation covariance. It is neither packet cadence nor six-state information.

**Failed pairwise relaxation.** Taking Delta(Qa) independently in norm costs
`min(2A_max,(J_max+Omega_max A_max)T_E)`. For every feasible theta_E<=pi/2,
both terms exceed `2g sin(theta_E/2)`: use theta_E<=Omega_max T_E,
`J_max+Omega_max A_max>g Omega_max` and `2A_max>sqrt(2)g`. Even before charging
bias and residual terms, this lower-margin estimate is nonpositive. This
invalidates the unsigned pairwise tactic, not the physical domain or the
signed whole-window target. Precision or independent boxes cannot repair it.

## A moving-history limitation of the physical span premise

The stronger intermediate claim that physical tilt span must produce changing
nominal attitude or nonzero innovations is false for some symbolic parameter
ranges. Let alpha=theta_E/2 in (0,pi/2), nu=2pi/T_E and
phi(t)=alpha sin(nu t). Take true body-to-world R(t)=Rx(phi(t)), p=v=a=0,
B=75e_x, physical gyro bias b_g=-omega and accelerometer bias
b_a=g(R(t)'e_z-e_z), with zero fast residuals. All measured values are exactly
level/north: f_m=-g e_z, gyro=0, mag=75e_x. The same deterministic shipping
execution is the quiet construction, with its actual references/gains/covariance.

Every T_E interval spans a full roll period and has gravity span 2alpha=theta_E;
it is not completely still. The physical biases are smooth and obey the
unchanged bounds under these sufficient symbolic conditions:

`g alpha<=B_a,s`, `g alpha nu<=D_a,s`,
`alpha nu<=min(B_g,s,Omega_max)`, `alpha nu^2<=D_g,s`.

Indeed ||b_a||=2g|sin(phi/2)|<=g alpha, ||dot b_a||<=g alpha nu,
||b_g||<=alpha nu and ||dot b_g||<=alpha nu^2. Translation, jerk and primitive
bounds hold exactly. The stationary actual-service lower bound applies at
each root multiplied by at least cos(alpha)^2 for the true heading/axial-bias
axis; require this remaining floor >1. These conditions have a nonempty
small-alpha, small-nu subfamily. They do not assign or qualify T_E/theta_E.

This exact conditional construction is not a counterexample to either
Delta_col or Delta_gyr: nominal force and field are nonparallel and nominal
gyro increment is zero. It invalidates only transfer through a mandatory
nominal tilt response, and shows why physical bias rates and signed identities
must remain even in the moving branch. It does not strengthen MAGNETIC SERVICE.

## Actual-S spline and chronological adjoint

For four actually applied S times t0<t1<t2<t3, put

`c_j=-6/product_(l!=j)(t_j-t_l)`,
`psi(t)=1/2 sum_j c_j (t-t_j)_+^2`.

The exterior jets psi, psi', psi'' vanish. Piecewise integration by parts
against the actual nominal v,p,S chain gives

`integral psi ahat dt=sum_j c_j r_S,j
 -sum_i (psi_i K_v,i-psi'_i K_p,i+psi''_i K_S,i)r_i-d_psi`.

At a spline atom the S correction uses the right psi'' and its pre-correction
S innovation. Same-time corrections before that S event use the left value;
later ones use the right value. The defect retains OU/source interpolation,
resets, frame changes, projection and arithmetic. No synthetic S observations
or independent innovation noise are introduced.

For `u_(i+1)=A_i u_i+K_i r_i+d_i`, arbitrary multipliers satisfy exactly

`sum W_i r_i=Z_N u_N-Z_0 u_0
 +sum (Z_i-Z_(i+1)A_i)u_i
 +sum (W_i-Z_(i+1)K_i)r_i-sum Z_(i+1)d_i`.

Use all needed mean coordinates; selecting BG/AW does not delete LIN/BA
coupling. Define chronological products

`Phi_(N,i+1)=A_(N-1)...A_(i+1)`, `C=[Phi_(N,i+1)K_i]_i`, `W=[W_i]_i`.

**Compatibility lemma.** The desired homogeneous adjoint exists iff
`ker C subset ker W`, equivalently one terminal multiplier satisfies `Z_N C=W`.
Backward propagation forces `Z_i=Z_N Phi_(N,i)`, proving necessity.
Conversely the map `Cx -> Wx` is well-defined on range C by the kernel
condition; extend it linearly and propagate backwards. No invertibility of
A_i is required. This retains actual rank-three gains and coordinate changes.
The exact rational checker returns either the multiplier or a witness
`Cv=0, Wv!=0`; a numerical rank tolerance never proves exact zero.

If Z_N=0, all Z_i=0. The first regular S atom has desired weight
`c_0(I-K_SS)`, nonsingular since `I-K_SS=R_eff(P_SS+R_eff)^-1` and R_eff>0.
The zero-terminal homogeneous specialization is therefore impossible on this
literal word. A forced adjoint moves forcing into the state residual sum;
it does not make that sum vanish.

`signed_temporal_diagnostic.py` tests the unrestricted terminal multiplier on
a quiet construction-to-release word, all 18 Euclidean mean coordinates and
four actual S events. Observer/control terminal parity is required. An exact
witness for exported rational coefficients concerns that carried word only,
not all-time magnetic service or a real-arithmetic trajectory enclosure.
The committed quiet word has knots (0.025,0.160,0.295,0.430) s relative to
its root and a seven-column witness with ||v||_infinity=1, Cv=0 exactly, and
|Wv|>754573637/5000000000000000. Its 80-digit value is approximately
3.0182945480072294e-7. This rejects exact homogeneous telescoping on that
exported word even with a free terminal multiplier. The actual quiet
innovations are zero; failure of the coefficient identity is not instability.

## Signed residual bounds and the remaining terms

Keep both compatibility residual sums. For example propagate Z backwards
homogeneously, leaving `E_i=W_i-Z_(i+1)K_i`, and substitute the literal acc,
mag and S identities into `sum E_i r_i`. Terms containing R_i ahat_i,
estimator biases, S means, the magnetic reference and literal defects remain.
Bound physical truth only after its temporal cancellations.

For matrix weights C_i and one physical bias history, exact finite summation
by parts gives

`sum_i C_i b_i=(sum_i C_i)b_0
 +sum_(j=0)^(N-2) (sum_(i>j) C_i)w_j`.

Its useful bound is

`B||sum C_i||+D sum_j ||sum_(i>j)C_i|| dt_j`.

This follows by substituting `b_i=b_0+sum_(j<i)w_j`. It charges only physical
rate when the signed total weight vanishes. It can be much smaller than
`B sum||C_i||`. It does not impose an estimator OU prior on physical truth.

For physical acceleration at cell starts, let C_i include the actual signed
weight and true rotation. The same-history velocity integral and jerk bound give

`sum_i h_i C_i a(t_i)
 =C_last v_N-C_0 v_0+sum_(i=1)^(N-1)(C_(i-1)-C_i)v_i+epsilon`,
`||epsilon||<=(J_max/2)sum_i ||C_i||h_i^2`.

Thus V_max multiplies endpoint norms plus signed weight variation; the
sampling error is charged afterwards. Nonuniform quadrature and the physical
displacement primitive retain their corresponding endpoint/difference terms.
The helper bounds require actual multiplier norm/variation ceilings. They do
not bound the remaining nominal terms or their correction action.

Zero-mean projection only cancels a constant physical bias in this Abel
identity. It generally destroys compatibility: `C=(1,2), W=(1,2)` is
compatible, but `W-mean(W)=(-1/2,1/2)` is nonzero on ker C. The former claim
`F_signed=F_defect+sum z_b w_g` without residual sums is withdrawn. Spline
boundary jets likewise do not solve the coupled mean adjoint. No estimator
BG or AW endpoint bound follows from either cancellation.

## Forced data adjoint on the same carried word

There is an exact alternative to assuming homogeneous compatibility. It
rewrites, rather than discards, both residual sums. In the default fixed
deployment frame use all 18 Euclidean means
`u=(bhat_g,vhat,phat,Shat,ahat,bhat_a)`. At a correction write

`r_i=y_i-Hbar_i u_i+epsilon_i`,
`u_(i+1)=A_i u_i+K_i r_i+d_i`.

Here Hbar is the **literal mean observation map**, not the full EKF Jacobian.
Its acc blocks are R_i at AW and I at BA; its S block is I at S; its magnetic
map is zero because actual attitude and reference remain in y. Consequently

`y_i^a=f_i^m+g R_i e_z`, `y_i^S=0`, `y_i^m=m_i^m-R_i B_i`.

The innovation arithmetic epsilon, additive mean arithmetic d, and separate
post-injection BA projection are retained. Predictions use the literal LIN
transport, BG identity and the mode-dependent BA predictor. Quaternion
corrections/resets affect subsequent actual R_i and K_i; no surrogate gains
or reset-free coefficient history is substituted. If additional configured
sensor-model terms or frame changes are present, their literal terms must
also be carried; they cannot be charged as independent sensor noise.

**Forced data identity.** For the existing desired innovation weights W_i,
set Z_N=0 and, in reverse chronological order,

`L_i=W_i+Z_(i+1)K_i`,
`Z_i=Z_(i+1)A_i-L_i Hbar_i`.

Then, with no rank or compatibility assumption,

`sum_i W_i r_i=Z_0 u_0+sum_i L_i(y_i+epsilon_i)+sum_i Z_(i+1)d_i`.

To prove it, put F_i=A_i-K_i Hbar_i. Then
`Z_i=Z_(i+1)F_i-W_i Hbar_i`; substitute the two literal recurrences and
sum the resulting adjacent state differences. The terminal term is zero
because Z_N=0; the root term generally is not. This is an inhomogeneous
adjoint. It does not repair the failed homogeneous equations or establish a
bound on its own action.

Using the same physical history, the acc and mag contributions become

`sum_acc L_i Q_i a_i`
`+sum_acc g L_i(R_i-Q_i)e_z`
`+sum_mag L_i(Q_i B_i^true-R_i B_i)`
`+sum_acc L_i(b_a,i+n_a,i)+sum_mag L_i(b_m,i+n_m,i)`.

Keep the two rotation/reference sums jointly signed. Apply the physical
velocity/primitive and bias-rate Abel identities to their respective
physical terms before taking norms. This identifies the missing quantities
precisely: root action, signed joint rotation/reference action, actual
multiplier norm/variation, and literal defects. Physical tilt span and
two-column magnetic information have not yet bounded them uniformly.
No independent nominal coefficient box is used. This lemma enters
`V_(j+1)<=rho V_j+c_d||d||^2` only through the unresolved signed separation,
six historical pivots and complete corrected-word matrix comparison.

**Gyro limitation of this identity.** In the default fixed-frame Euclidean
recurrence, `A_i E_BG=E_BG` and `Hbar_i E_BG=0`, hence
`F_i E_BG=E_BG`. Induction gives `Z_i E_BG=0` for the zero-terminal forced
innovation functional, regardless of the BG rows of the actual gains.
Thus this identity removes its BG root coefficient but supplies no restoring
feedback for the BG mean itself. This unit column is a property of the
frozen-attitude algebra, not an unstable mode of the coupled nonlinear
estimator. The shipping projection independently supplies the BG level bound.
The full signed gyro margin still needs actual attitude recurrence and
magnetic/tilt constraints; it cannot follow from this Euclidean matrix alone.

**Moving carried diagnostic and one refinement.** The observer runs the
existing smooth source from construction: roll `0.02 sin(t/2)` and vertical
displacement `0.4 sin(0.6t)`. Its physical bounds are satisfied; every window
of length at least `4 pi` has gravity span `0.04`, so its attitude premise
is conditional on `T_E>=4 pi`, `0<theta_E<=0.04`. This is not a numerical
qualification of the symbolic deployment parameters. All-time actual
MAGNETIC SERVICE remains unverified for this finite diagnostic.

The 225--289-s carried tail uses four actual S knots
`(0.065,21.365,42.665,63.96)` s relative to the root. The balance includes
the inherited mean at 225.06 s through the last correction at 288.96 s;
the 225-s origin is used only to label the spline knots. There is no reset.
The balance includes
12780 predictions, 12780 acc, 1598 mag and 892 S corrections, and 15270
post-injection projection observations. All 18 means have consecutive
operation-boundary continuity. The observer and untapped control have
identical terminal means, quaternion, covariance, stages and magnetic count.
The actual reference is checked constant on this particular word.
The declared physical gravity 9.80665 and its shipping float representation
are kept separately; their difference is an explicit model/arithmetic term.
At 80 digits the forced identity residual norm is approximately
`2.84874602317e-82`; this is a finite decimal diagnostic, not an interval
enclosure or a source-uniform certificate. The root multiplier is not zero:
its block Frobenius norms for BG,V,P,S,AW,BA are approximately
`(0,.01746150,.003428685,.001083323,.00008579003,1.117579)`.
Only its BG cancellation is the structural identity proved above.

The recorded reference gives projected gravity `g||Pi_B e_z||=8.77133455729`.
This is a diagnostic necessary budget for the normalized actual-S average
under nominal collinearity; the spline defect and all other supplies still
have to be charged. It is not a source-uniform nominal field floor. The
rotation triangle bound costs `39.22565052489`, giving margin
`-30.45431596760` before all other charges. This invalidates that norm
relaxation, not the actual signed rotation term (whose norm is about .00116)
or the original separation target. No interval enclosure of it is attempted.

The cellwise signed-velocity/jerk bound costs
`239.57211610857 + .39841221469 = 239.97052832325`. One refinement combines
signed coefficients C_i inside each actual S interval [a,b] first. Put
H=b-a and D=sum C_i. The exact physical identity and remainder are

`sum C_i a(t_i)=(D/H)(v(b)-v(a))+epsilon`,
`||epsilon||<=J sum ||C_i||[(t_i-a)^2+(b-t_i)^2]/(2H)`.

Proof: compare each a(t_i) with `H^-1 integral_a^b a(s) ds` and integrate
`||a(t_i)-a(s)||<=J|t_i-s|`. The velocity endpoints are shared across
intervals, so sum by parts again with the signed D/H before taking norms.
For C=(1,-1), sample times (0,1), a(t)=t, the total D vanishes but the
remainder is exactly -1, attaining the bound; it cannot be discarded.

On this word the refined velocity charge is `23.98805974096`, jerk charge
`3.82238262046`, total `27.81044236142`, still `19.03910780414` above the
recorded threshold before root, sensor, bias or spline defects. The signed
physical acceleration action itself has norm about .005, but that value
cannot be used as an upper bound for other histories. After this second
failure, stop the variation-norm relaxation. The architecture review retains
the same construction-to-H18-to-A21 path: the missing step is a joint bound
on physical integral action, rotation/reference action and chronological
gains, together with the attitude recurrence for BG level. Neither a nominal
box nor a new proof track supplies it. Finite diagnostics alone cannot test
the universal implication or justify strengthening MAGNETIC SERVICE.

## Construction and gyro complete turns

The MEKF gyro-bias mean begins at zero; Live handoff imports no proxy bias,
and pre-prediction magnetic updates have zero BG gain. The first increment is
at most `.006*(.6108652381980153+.02+.02)=.0039051914291880918` rad.

Rodrigues integrated bias transport has transverse singular values
`h |sinc(theta/2)|`, vanishing at nonzero complete turns. Subsequently
`omega_hat=omega+b_g+n_g-bhat_g`; a complete turn at any allowed step requires

`||bhat_g||>=2pi/h_max-(Omega_max+B_g,s+B_g,f)>1046.5466859583 rad/s`.

The shipping residual gyro-bias mean now has the implemented Euclidean bound
`||bhat_g||<=.5 rad/s`, independently of the physical `.02 rad/s` bias bound.
This excludes that barrier from construction and every subsequent qualified
prediction. The complete rate envelope gives `h||omega_hat||<=.0069051914291880918`
rad; including the normalized quaternion polynomial remainder gives `.007` rad.
Both source branches have a transport singular floor `.003999991833333333 s`.
See `ou-gyro-bias-projection.md` for the exact rational certificate and the
limitation to qualified 4--6 ms steps, rather than arbitrary API-admitted stalls.

The historical `h=.005, bhat_g=-400pi e_z` example violates the implemented
invariant. This does not close the full signed `Delta_gyr(W)` or a historical
six-pivot bound: actual chronological observations, resets and projection
defects still require their joint enclosure. The projection changes only
the mean. Its ideal component sector has gap `R_g-B_g,s=.48 rad/s`; covariance
cross terms and inward floating-point rounding remain separate signed supplies.
The forced balance observes both gyro and accelerometer mean projections.

## Historical reader: the corrected conditional implication

The previous assertion that Delta_col/Delta_gyr exclude every six-column
rank-loss mechanism, with pivots >=min(Delta_col,Delta_gyr), lacked a proof
and is withdrawn. Coefficient compactness also needs a derived nominal-history
bound. The physical 3-D vector Gramian is not the historical observation array.

The valid conditional algebra starts with all SIX unsquared Gram--Schmidt
residual norms of the existing largest-residual selected minor O_I bounded
below by p>0. The selector compares their squares, which has the same order.
If C>=1
bounds operator norms of that minor, transports and H_i, and H>=1 bounds
both the terminal selector and T_h, determinant/singular-value comparison gives

`||O_I^-1||<=C^5/p^6`, `||L||<=ell:=H C^5/p^6`, `L O=T_h`.

For at most N operations the backward recursion, including its nonzero
terminal selector and all `Y<-Y-L_i H_i` updates, gives

`||Y_i||<=Ymax:=C^N(H+N ell C)`.

If E bounds complete correlated process/measurement factors and U bounds the
inherited nuisance covariance, exact AG root cancellation gives

`B_W<=[N E^2(Ymax^2+ell^2)+U Ymax^2]I_6`.

This coarse finiteness lemma does not replace the full matrix process
comparison by a scalar contraction bound. No uniform p,C,H,E,U,N is asserted
for the actual history family. The margins must first be proved and linked
to the six pivots. Only then can the established implication
`B_* -> AG upper -> J_AG>0 -> full upper -> rho_0<1` be instantiated.
No new contraction enclosure is attempted without feasible full-matrix input.

## Controlling unresolved implication

The retained physical reserve is exactly
`27049050188592/625000000000000000=4.32784803017472e-5`, under the constant
physical field hypotheses of its sampling lemma. Its nominal transfer is
not bounded by physical tilt span alone. The unresolved quantity is that
reserve minus compatible signed transfer, compatibility residuals and
literal defects over every carried history. Neither Delta_col nor Delta_gyr
has a positive certified uniform bound; neither has been rigorously falsified
under all current assumptions.

The smooth diagonal-wave R=I diagnostic is outside the strengthened moving
domain. Quiet water and quiet-bias ambiguity remain inside; retain the literal
bias-projection sector rather than universal V<=36 entry. The nonlinear
radius, capture, H18/release, magnetic-service qualification, every-prefix
retention and float32 bounds remain open. The unsigned cumulative-energy
ellipsoid's approximately 8e5 deficit is a DEAD_END and is not revisited.
