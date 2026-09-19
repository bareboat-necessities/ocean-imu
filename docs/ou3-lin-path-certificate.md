# LIN endpoint-path factor certificate

## Role in the single proof

This is the constructive LIN covariance lemma within the AG/LIN/BA factor
metric. It does not create a second theorem path. The long-word neutral
information floor remains `mu_N >= 2.04e-3`; covariance normalization, finite
error, arithmetic and capture remain separate fail-closed obligations.

The executed rational comparison gives `ell_LIN >= 1.196007314542406680e-6`
after the complete correction window, conditional on the stated shipping
noise-floor binding and real-arithmetic covariance equations. It covers every
piecewise-constant `tau(t) in [0.02,12]`, every mesh step in `[0.004,0.006]`, and
every correction subset. No frozen-tau sampling or one-step scalar q_min is
used. Float32 transfer is NOT claimed by this result.

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
Taking the trace of each PSD endpoint-energy matrix gives the certified
precision ceiling and its reciprocal factor floor.

The deployed integral standard deviation has the lateral factor 0.50, not
0.72 on both axes. Hence the comparison uses `0.15*0.50=0.075`, not 0.108.
The initial unbound calculation at 0.108 was conditional, not a deployment
certificate. Noise-floor binding is kept open until the complete configuration
and parameter-smoothing paths are audited.

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

The resulting posterior precision ceiling is at most
`6.990887744796e11`, with covariance floor greater than `1.4304334964e-12`
in the retained `(5.5,8.1,1100,4)` LIN coordinates. All certificate arithmetic
is rational; decimal endpoints and the square-root lower export are outward.

## Remaining falsifiable transfer

Bind all source polynomial coefficients and deployed noise floors, including
smoothing, hard events and anisotropy. Then enclose the literal float32
covariance perturbation in this factor metric. A positive real-arithmetic
factor is not permission to set constructive_root_covariance_floor or
constructive_full_A21_mu_rho_enclosure true. No numerical cross-block ceiling,
complete rho0, nonlinear radius, startup capture, or practical radius is
claimed by this module.


## Promotion audit

The interleaved-correction audit distinguishes this variational certificate from
the retired generic endpoint helper.  The latter attempted to infer a posterior
floor from only an accumulated endpoint process floor plus event information
ceilings; that implication is false.  Here the comparison is instead an
endpoint minimum-action problem: each trial path pays its process action and
the measurement action at the literal event times.  The measurement mesh bound
therefore applies before minimization and directly upper-bounds endpoint
precision.  The two-dimensional counterexample to the generic helper is
admitted by this formulation and does not contradict the action inequality.

Shipping bindings were checked against the authoritative implementation:
`tau` is clamped to [0.02,12] s; the a_w stationary standard deviation is
floored at 0.05 m/s^2; the deployed accelerometer standard deviation is 0.2
m/s^2 before only nonnegative vibration/noise inflation, so the certificate's
0.05 lower bound is conservative; the S base standard deviation is clamped at
0.15 m*s and the smallest deployed axis factor is 0.50, giving the used 0.075
floor.  The default periodic a_w synchronization is a queued PSD covariance
inflation inside prediction, so omitting it can only reduce the comparison
covariance.  The real-arithmetic shipping factor is therefore promoted as
`ell_LIN >= 1.196007314542406680e-6`.

Combining this with `ell_AG=4.999974644e-4`,
`ell_BA=5.618273739e-4`, root `gamma=1`, and
`mu_N>=2.04e-3` gives the conservative scalar normalization
`mu_cov>=2.9180843327e-15` and
`rho0<=0.9999999999999971`.  This closes positivity of the constructive
real-arithmetic linear certificate, but the margin is extremely small.  The
next falsifiable step is the explicit nonlinear remainder and float32 supply:
they must fit inside this margin, or a sharper factor-coordinate normalization
must be derived without weakening assumptions.


## Matched factor-coordinate sharpening

The previous normalization used the retained physical-envelope coordinate scales
(5.5,8.1,1100,4) and then collapsed the complete LIN factor to its smallest
singular value.  That collapse is coordinate dependent.  A fixed diagonal
congruence is therefore allowed provided *both* covariance and information are
recertified in the same coordinates; changing the estimator, motion contract,
or measured data is not allowed.

The declared proof coordinates are now tested at
`(v,p,S,a_w)=(2.4,18,132,4)`.  These are proof coordinates only, not tighter
physical bounds.  The exact-rational path-action certificate is rerun with that
map, and the analytic three-S-row information certificate is independently
rerun with the same map.  It gives `mu_trans > 2.1e-3`, while the LIN factor
is about `8.91e-6`.  The resulting conservative neutral normalization is

`mu_cov >= 1.66e-13`,
`rho0 <= 0.9999999999998335`.

This is more than fifty times the prior scalar-coordinate margin, without
sampling trajectories or fitting rho.  It is nevertheless far too small to
assume that nonlinear and float32 supplies will fit.  The remaining avoidable
loss is the scalar collapse itself.  The next sharpening must retain the full
LIN endpoint action/factor matrix and the compatible translation information
matrix and certify `lambda_min(L^T J L)` directly.  No physical assumption or
quality gate is to be tightened if that matrix certificate is still
insufficient.
