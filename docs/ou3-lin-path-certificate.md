# LIN endpoint-path matrix certificate

## Role in the single proof

The exact-rational LIN comparison supplies a covariance factor for the
complete-word energy inequality. It certifies the full 4x4 action matrix A and
its inverse in coordinates (v,p,S,a_w) scaled by (2.4,18,132,4). The factor is
not reduced to a smallest LDL pivot or an independently normalized scalar.
The companion joint-root certificate combines this full Loewner comparison
with fresh AG/BA process injection at the same post-prediction root.

The result covers arbitrary piecewise-constant tau in [0.02,12], mesh steps in
[0.004,0.006], all regular interleaved corrections, and default PSD covariance
sync. The implementation comparison is real arithmetic; float32 transfer and
full-state contraction remain open.

## Variational factor construction

Write a scalar LIN path as `(v,p,S,a)=(S'',S',S,S''')`. A degree-seven Hermite
polynomial has four zero initial jets and any prescribed final jets. This is a
comparison path for a covariance inequality; neither the actual estimator nor
its physical history is reset. At any endpoint after 16 seconds, choose the
sample just before its 16-second look-back, so the comparison duration belongs
to `[16,16.006]` seconds. Discarding the actual PSD prior only reduces the
comparison covariance.

For `lambda=1/tau`, the path's control is `a'+lambda a`. With
`q(t)>=2 sigma_min^2 lambda(t)`, its process action obeys

`integral (a'+lambda a)^2/q <= [12 integral (a')^2 + 50 integral a^2 + a(T)^2] / (2 sigma_min^2)`.

The cross term telescopes without differentiating lambda. Thus arbitrary tuner
changes are admitted. A trial path upper-bounds the minimum endpoint action,
which is the inverse controllability Gramian quadratic form.

## Literal correction survival

For an accelerometer row partitioned into AG/LIN/BA,
`H' R^-1 H <= 3 diag(H_AG' R^-1 H_AG,H_LIN' R^-1 H_LIN,H_BA' R^-1 H_BA)`.
This introduces more informative block-separated observations, so Riccati
monotonicity gives a lower covariance comparison. The LIN accelerometer part
observes only `a` through an orthogonal rotation. Integral observations select
`S`. Applied magnetic observations affect the AG comparison only; no rejected
or unavailable event is counted as magnetic information.

For a sample mesh with separation at least delta and any path with f(0)=0,

`sum_i |f(t_i)|^2 <= (2/delta) integral |f|^2 + 2 delta integral |f'|^2`.

Proof: compare f(t_i) to f(t) on disjoint backward intervals of length delta,
use Cauchy-Schwarz and `(x+y)^2<=2x^2+2y^2`, and integrate. Apply this to both
`a` and `S`. It covers an integral update at every IMU sample, hence every
literal sparser scheduler history, and includes accelerometer corrections at
every sample. Add the resulting measurement energies to the process action.
The variational posterior precision is no greater than this trial action.
Retain the complete 4x4 endpoint-energy matrix A. The variational
comparison gives E_LIN^T P^-1 E_LIN <= A and therefore the full-state
Loewner lower bound P >= E_LIN A^-1 E_LIN^T. This conclusion includes
carried cross covariance; it is stronger than a marginal P_LL bound.

## Literal small-argument polynomial defects

The shipping OU transition and covariance use finite polynomials for
`x=dt/tau<0.01`. These are not silently identified with an exact exponential.
In natural step coordinates `D_h=diag(h,h^2,h^3,1)`, write
`Q_exact=sigma^2*x*D_h B(x) D_h'`. For integration orders m,n in `(1,2,3,0)`,

`B_mn(x)=2 sum_r (-x)^r sum_{j=0}^r [1/((m+j)!(n+r-j)!(m+n+r+1))]`.

The retained degree is read from the literal h^9 covariance formulas. A
positive-series geometric tail majorant and the zero-drift controllability
inverse give the rational relative defect ceiling
`epsilon_Q < 1.944803e-5`; thus `Q_poly >= (1-epsilon_Q) Q_exact`.
The transition defects obey
`|Delta F_pa|<=h^2*x^3/120` and `|Delta F_Sa|<=h^3*x^3/720`.
Their endpoint action is bounded using the same natural-step inverse and added
to the process action via Young's inequality. This is a relative factor
comparison, not a tiny raw one-step eigenvalue replacement.

## Matrix enclosure and exact verification

For each derivative order, the Hermite Gram matrix is evaluated at 16 s.
Every entry is a constant times an integer power of the duration, so endpoint
values at 16 and 16.006 s bound its variation. Adding the maximum absolute
row sum of that variation times I gives a uniform Loewner upper matrix.
Process, source-defect and measurement matrices are then combined before
inversion. `lin_matrix_certificate.py` verifies exact positive LDL pivots and
exports A and A^-1 in `lin-matrix-certificate.json`.

A 70-digit non-promoting diagnostic preceded the rational enclosure. Comparing
the retained action matrix with raw neutral S rows at actual event boxes
[0,.156], [8,8.156], [16,16.156] gives a lower normalization of
3.5959862602014e-11. The adjugate/determinant comparison is rational and covers
the full boxes. This is not the transported closed-loop loss matrix and does
not yield a full A21 rho. No independent heading information is inferred from
the restricted magnetic-service Gramian.

## Bound shipping profile

The deployed default SpectralMSE law has unit cadence renormalization and
S_factor=1. The a_w stationary standard deviation is floored at .05; S base
standard deviation is clamped at .15 with smallest axis factor .50, giving
.075. The deployed accelerometer nominal standard deviation is .12 before
nonnegative inflation, so the .05 comparison floor is conservative. Periodic
a_w synchronization queues a PSD increment inside prediction; omitting that
increment reduces the comparison covariance. Alternate policies and arbitrary
user setters are outside this binding.

## Joint post-prediction lower bound

At the same regular A21 root, fresh process noise gives
`P >= X=diag(q_AG I_6,0,q_BA I_3)`, while the LIN comparison gives
`P >= Y=diag(0,A^-1 tensor I_3,0)`. Hence `P >= (X+Y)/2 > 0`.
The exact combined AG and BA diagonal floors are respectively 4.9991e-10 and
at least 4.99999e-10; the LIN block stays A^-1/2. All 21 coordinates are covered.
The comparisons retain cross covariance and do not add unrelated marginal floors.
`root_covariance_certificate.py` binds the fresh AG/BA source constants and
verifies the combined factors by exact rational arithmetic.

This closes the constructive real-arithmetic lower covariance at the specified
roots. It does not establish uniform upper covariance, every-prefix retention,
full-state information/loss, a nonlinear radius, capture, or float32 stability.
The full theorem and constructive_full_A21_mu_rho_enclosure remain open.
