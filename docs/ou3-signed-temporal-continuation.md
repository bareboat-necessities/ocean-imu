# OU-III signed temporal continuation

This note records the active source-uniform construction for the single shipping-faithful stability proof. It does not change estimator behavior, physical assumptions, quality gates, or theorem architecture.

## Signed temporal identity

Choose four strictly increasing times t0<t1<t2<t3 at actually applied integral corrections and define

c_j = -6 / product_{l != j}(t_j-t_l),
psi(t) = 1/2 sum_j c_j (t-t_j)_+^2.

Then D^3 psi = sum_j c_j delta(t-t_j), with psi, psi', psi'' vanishing outside [t0,t3]. Piecewise integration by parts against the literal nominal chain, retaining every correction jump, gives

int psi a_hat dt
 = sum_j c_j r_S,j
 - sum_k (psi_k K_v,k - psi'_k K_p,k + psi''_k K_S,k) r_k
 - d_psi.

Here r_k is the chronological innovation of the actual applied accelerometer, magnetic, or integral correction, and d_psi retains source-transition, hard-event, reset/interpolation, and arithmetic defects. The identity does not manufacture pseudo-observations between callbacks and does not drop the actual v,p,S jumps.

For the affine mean recursion u_{i+1}=A_i u_i+K_i r_i+d_i, summation by parts gives

sum_i W_i r_i =
 Z_N u_N-Z_0 u_0
 +sum_i (Z_i-Z_{i+1}A_i)u_i
 +sum_i (W_i-Z_{i+1}K_i)r_i
 -sum_i Z_{i+1}d_i.

Exact endpoint elimination therefore requires both Z_i=Z_{i+1}A_i and W_i=Z_{i+1}K_i. A zero-terminal homogeneous adjoint cannot represent the first nonzero integral atom: Z_N=0 forces all Z_i and W_i to zero, while c_0(I-K_SS) is nonsingular at a regular applied S correction. The proof must therefore use an observation-forced adjoint and retain any remaining signed innovation functional explicitly. This is a failed multiplier specialization, not a physical counterexample.

## Physical temporal enclosure

Under the unchanged regular timing and MARINE MOTION bounds, the four-event weight yields the exact rational sampled-acceleration mean ceiling

1050297/4840000 = 0.217003512396694... m/s^2.

Combining the same-history acceleration balance with the existing magnetic residual limits gives a physical joint-vector information floor

gamma_phys > 2.432784801508736e-4,

with margin >4.32784801508736e-5 above 1/5000. This is physical measured-vector information only; it is not promoted to nominal AG information or historical-reader rank.

## Projection sector

For the literal radial accelerometer-bias projection, d_b=b_hat_corr-Proj(b_hat_corr), physical ||b_a||<=B_a and projection radius R_b>B_a imply

(e_b_corr)^T d_b + ||d_b||^2 + (R_b-B_a)||d_b|| <= 0.

For H=P^{-1}, choose lambda>0 and D_b=2 lambda I-E_b^T H E_b positive definite. With h_b=(E_b^T H-lambda E_b^T)e, completing the square gives

V(e+E_b d_b)-V(e)
 <= (h_b-y)^T D_b^{-1}(h_b-y)

for every ||y||<=lambda(R_b-B_a). Thus ||h_b||<=lambda(R_b-B_a) is an explicit shaped projection-nonexpansion condition retaining all covariance cross terms. Outside it the projection defect remains a signed supply. Universal entry into V<=36 is neither required nor asserted.

The unchanged bias constants give R_b-B_a > 0.1748333950160459 m/s^2. A regular-A21 pre-projection BA precision ceiling 1000003000 I3 is available after charging the binary32 driving coefficient; this is an operation bound, not an invariant-region certificate.

## Positive temporal margins imply a common historical action ceiling

The next implication is now closed analytically. Normalize the six AG root
coordinates as in the temporal margins and put

delta = min(inf_W Delta_col(W), inf_W Delta_gyr(W)) > 0.

On the separated same-history coefficient family, use the existing
largest-residual factor pivot rule on the complete six-column raw observation
array O(W). The two temporal exclusions rule out the only sustained AG rank-loss
mechanisms left by the literal attitude/bias transport: force/field
collinearity and complete-turn gyro aliasing. Therefore each W has a full-rank
six-row pivot chart. Strict delta separation and compact shipping
coefficient/factor bounds make the residual pivot functions continuous and
bounded away from zero on a finite chart cover.

More quantitatively, if C bounds the normalized observation coefficients and
each of the six successive residual pivots is at least delta, the selected
minor O_I obeys, by the determinant product and adjugate/Hadamard bound,

||O_I^{-1}||_2 <= C^5 / delta^6.

For H bounding the terminal AG map, the exact reader L=T_h O_I^{-1} therefore
satisfies

||L||_2 <= H C^5 / delta^6.

The historical backward recursion contains a fixed finite number N of
prediction/reset transports and rank-at-most-three observation/process
factors. If C also bounds each chronological transport coefficient, E bounds
the corresponding noise/process factors, and U bounds the inherited nuisance
covariance, then the exact matrix action satisfies the explicit conservative
ceiling

B_W <= B_* I_6,

B_* = N (H C^5 delta^-6 E C^N)^2
      + (H C^5 delta^-6 C^N)^2 U < infinity.

This keeps the coefficient-dependent rank-three reader/factor structure and
does not invoke a scalar information lifting or Riccati box. The scalar display
is only an operator-norm ceiling on the already full matrix backward action;
the implementation retains the exact matrix action for the later J/rho
enclosure.

Hence

inf_W Delta_col(W)>0 and inf_W Delta_gyr(W)>0
  ==> L(W)O(W)=T_h(W), B_W <= B_* I_6 < infinity.

The implication is theorem-level algebra. Its premises are not yet promoted:
the source-uniform outward enclosure of the two signed temporal margins remains
the controlling numerical/analytic task.

## Executed source-uniform margin attempt

The enclosure attempt is now executed rather than deferred.  The physical
collinearity reserve available to the signed balance is the exact rational

27049050188592/625000000000000000
 = 4.32784803017472e-5.

To promote Delta_col this positive reserve must exceed the operator norm of
the observation-forced signed innovation transfer plus the literal
reset/transition/hard/arithmetic transfer on the same carried window.  The
unchanged assumptions and currently proved covariance statements do not yet
supply a source-uniform ceiling for that forced-adjoint functional.  Assigning
one from the finite construction trace, NIS energy, or independent gain boxes
would revive a forbidden relaxation.  The attempted inequality therefore
stops exactly at

physical reserve - signed forced-adjoint/reset transfer > 0,

with the second term unbounded by the present theorem premises.

The gyro balance reaches the same structural point. MAGNETIC SERVICE supplies
normalized recurring root information, and physical gyro bias/rate increments
are bounded, but conversion to a lower bound on the two transverse singular
values of the literal integrated bias transport requires a source-uniform
bound on the signed corrections that move b_hat_g and on reset transport.
That bound is likewise not presently derivable from the retained lemmas
without using unsigned innovation energy or an independent nominal box.

Thus this PR has executed the requested enclosure attempt and identified the
first failed inequality. No negative Delta margin is claimed: the missing
quantity is a rigorous ceiling, not a constructed adverse history. No new
physical assumption is introduced.

## Literal forced-adjoint source bound

The observation-forced correction functional can be reduced further without
an innovation-energy estimate.  For the literal affine mean recursion and a
compatible adjoint,

sum_i W_i r_i =
 Z_N u_N - Z_0 u_0 - sum_i Z_(i+1) d_i.

Hence, if zeta bounds the relevant endpoint multipliers, U bounds the
corresponding nominal endpoint mean and D bounds the accumulated literal
affine defects,

|sum_i W_i r_i| <= 2 zeta U + sum_i ||Z_(i+1)|| ||d_i||.

This is the requested source-uniform form: innovations and gains disappear
through the signed chronological identity rather than being independently
bounded.

Applying it to the shipping source exposes a sharper obstruction.  The
accelerometer-bias estimate has the literal global projection bound
||b_hat_a||<=0.4.  In contrast, the gyro-bias estimate has identity mean
prediction and no projection/saturation, and the latent AW mean has OU
prediction/corrections but no mean saturation.  The physical bounds
||b_g||<=0.02 and ||a||<=8.8 do not bound those two estimator means before the
same-history error/capture result being proved.  Covariance/tuner clamps are
not mean-state clamps.

Therefore the exact telescoping derivation removes the previously missing
innovation-functional ceiling, but replaces it by BG/AW endpoint terms for
which the current theorem premises contain no absolute source-uniform bound.
Using the finite 600-s mean audit (BG<1), or imposing a nominal AW box, would
be circular/nonuniform and is not promoted.

This is an analytical endpoint obstruction, not a numerical-conditioning
failure: under the unchanged theorem premises the present forced-adjoint
choice cannot yield a finite numeric source-uniform correction/reset ceiling.
A successful next construction must choose multipliers/physical balances whose
BG/AW endpoint coefficients vanish, or cancel those endpoint means against
the gyro and integral physical identities. Adding an estimator-mean bound as
a new physical assumption is not justified.

## Endpoint-annihilating multipliers

The uncontrolled nominal endpoint means can be removed analytically rather
than bounded.

For the translation/AW chain, repeated integration by parts leaves AW boundary
coefficients proportional to psi, psi' and psi''.  The four-actual-S-event
quadratic spline already satisfies

psi=psi'=psi''=0

at both ends of its support.  Consequently the signed physical v,p,S identity
has no absolute a_hat_w endpoint term.  AW enters only through the chronological
correction/source defects already retained by the identity.

For the attitude/gyro-bias pair, let z_theta be the left temporal multiplier
transported by the literal attitude adjoint. Introduce a companion multiplier

dot z_b = -z_theta.

Set z_b(t0)=0.  The terminal absolute gyro-bias coefficient also vanishes
exactly iff

integral_[t0,t1] z_theta(t) dt = 0.

Thus the needed gyro construction is the zero-mean subspace of transported
attitude multipliers.  On a discrete carried word, if q_k is the exact cell
integral of z_theta, replace a candidate sequence by

q_k^0 = q_k - (1/N) sum_j q_j.

Then sum q_k^0=0 and the companion recursion
z_b,k+1=z_b,k-q_k^0 starts and ends at zero.  Applying the signed physical
gyro-bias recurrence leaves only

sum_k z_b,k w_g,k,

with ||w_g,k||<=D_g dt_k, plus the literal correction/reset/arithmetic defects.
There is no absolute b_hat_g endpoint term.

After both endpoint cancellations the source functional has the form

F_signed = F_defect + sum_k z_b,k w_g,k,

and therefore

||F_signed|| <=
 sum_k ||Z_(k+1)|| ||d_k||
 + D_g sum_k ||z_b,k|| dt_k.

This is source-uniform once the transported multiplier norms and literal
defect magnitudes are enclosed; neither a nominal AW box nor an absolute BG
estimate bound is required.

The remaining issue is no longer endpoint boundedness.  The zero-mean
projection must be composed with the actual magnetic/accelerometer observation
forcing and shown to retain a strictly positive collinearity and gyro
separation margin after the multiplier/defect norm cost.  That is the next
quantitative inequality.

## Controlling unresolved implication

The first unresolved mathematical implication is source-uniform control, on every admitted carried same-history window, of the remaining signed acc/mag/S innovation functional together with construction-linked gyro-bias evolution and literal transition/reset/hard/arithmetic defects, with enough strict margin to prove

inf_W Delta_col(W)>0,  inf_W Delta_gyr(W)>0.

The implication from those two margins to L(W)O(W)=T_h(W) and a common finite B_* is now proved above. What remains unproved is the source-uniform positivity of the two premises themselves; B_* therefore remains uninstantiated and is not promoted as a shipping certificate. Consequently J_AG, the full 21-state covariance upper bound, rho_0<1, the complete nonlinear retained region, capture/release entry, every-prefix retention, recurring source-qualified magnetic service, and float32 totality remain open.

Do not replace this obligation with independent nominal boxes, unsigned innovation energy, restricted-information lifting, Riccati subdivision, scalar/Gershgorin reductions, or finite-history diagnostics.
