# Full-state loss factors and nuisance elimination

This construction enters the existing tail inequality through
`D_W >= delta P_0^-1`, delta>0. It changes proof computation only. The dense
exact-rational checker remains an independent reference on adversarial words.
The source-uniform full-state margin remains open.

## Measurement factor and covariance array

Write `P=C C^T`, `P_0=C_0 C_0^T`, and carry `T=M C_0`. For an applied
observation, let `S=L_s L_s^T=H P H^T+R_eff` and define

`B_m=L_s^-1 H T`, `T^+=T-K(H T)`.

Then `B_m^T B_m` is the root-whitened loss. The triangular solve has at most
three rows. Retain B_m until streaming QR compression; do not construct its
21x21 normal matrix at each event. For covariance, orthogonally triangularize

`[[R_eff^(1/2), H C], [0, C]]`.

The resulting lower array has blocks `[[L_s,0],[K L_s,C^+]]`. Its bottom
Schur factor is the exact Joseph posterior in real arithmetic. Cross-state
terms are preserved, including the nuisance directions in H T.

## Prediction factor without precision subtraction

For a supplied rectangular process factor U, `Q=U U^T`, set `W=[F C,U]`.
Use the complete QR factorization

`W^T = [O_1,O_2] [L^T;0]`.

The covariance factor after prediction is L. Let J select the first n
coordinates of the n+r dimensional orthogonal array. Then

`A=L^-1 F C=O_1^T J`, `A^T A+(O_2^T J)^T(O_2^T J)=I`.

Consequently the prediction loss factor is

`B_p=O_2^T J C^-1 T`.

This factors the decrement by an orthogonal complement, with no subtraction
of nearly equal precision matrices. It admits singular Q; F need not be
invertible if the predicted covariance is positive definite. A nonsingular
reset replaces C by a triangular factor of G C and T by G T, adding no loss.

The engine appends B_m and B_p by QR, retaining at most n=21 rows in R_W.
At the endpoint,

`Z=L_W^-1 T_W`, `Z^T Z+R_W^T R_W=I`.

Thus the exact implication is `sigma_min(R_W)^2>=delta => rho0<=1-delta`.
Products used to inspect this identity are diagnostics, not the accumulation
representation. Three sensor channels do not make the complete covariance
rank three. The magnetometer skew block has rank two; integrated OU process
blocks have four controllable coordinates per axis.

## Nuisance information

Partition a loss factor by state columns as `[B_h,B_n]`. Information remaining
after nuisance adjustment is

`min_n ||B_h h+B_n n||^2 = ||(I-Proj_range(B_n)) B_h h||^2`.

When B_n has full column rank, this is the ordinary Schur complement. In the
singular case, choose an exact independent-column basis N of its range, solve
`(N^T N) X=N^T B_h`, and retain `B_h-N X`. Exact orthogonality to *all*
nuisance columns certifies the projection. No singular inverse or numerical
rank cutoff is used in the rational certificate.

The audit factor `[[1,0,1,1],[0,1,0,0]]` has a singular nuisance Gramian and
reduced heading information diag(0,1). Adding row `[0,0,1,1]` changes that
reduction to diag(1/2,1), but leaves a nuisance nullspace: the full factor
still has rank three with four state columns. A positive heading Schur
complement is therefore not itself a full-state contraction certificate.

For a square rational factor B and candidate right inverse X, the verified
residual bound is

`sigma_min(B_true) >= (1-||B X-I||)/||X|| - eta`,

provided `||B X-I||<1` and `||B_true-B||<=eta`. The implementation uses exact
rational products, Frobenius upper bounds and directed square-root bounds.
Proving eta for the *whole word* is mandatory. Setting eta=0 certifies only
the stored rational matrix, not the preceding floating-point QR computation.

## Executed feasibility and source parity

The exporter instantiates the unchanged shipping core with double arithmetic,
constant quiet-water data, zero physical biases, field (20,0,40) microtesla,
dt=.005, tau=1.8, sigma_aw=.05, S standard deviations (.108,.075,.15),
acceleration std .12 and magnetic std .25. After 3200 carried warmup samples,
it exports the next 3200 samples, without state/covariance/scheduler reseeding.
Every magnetic/accelerometer update must report acceptance. The source state
and attitude stay exactly zero/identity; all actual injection resets are
therefore identity. This is a configured core diagnostic, not a wrapper
startup/capture or adaptive-tuner qualification.

The 16-s word contains 3200 predictions and 4907 corrections. The measured
full-state values are delta≈.1672905911 and rho≈.8327094089. The identity
residual is 1.37e-12; maximum relative source-double gain/covariance mismatch
is below 5e-14. A 70-digit independent dense calculation on 60 prefix events
gives rho=.9999887113430517, agreeing with the factor computation. This check
preceded the rational residual enclosure of the stored factor matrix.

The 16 disjoint one-second magnetic restrictions have measured minimum above
13.5 in the declared heading/gyro scales. This does not prove every sliding
window or an all-time service certificate. Full nuisance elimination is
performed before any heading diagnostic is reported. No number here is
promoted to a source-uniform theorem rho or float32 bound.

Reproduce from the repository root (set EIGEN_DIR to an installed Eigen root):

```sh
g++ -std=c++17 -O2 -I src -I "$EIGEN_DIR" tools/stability/ou3_theorem/export_quiet_word.cpp -o /tmp/ou3-export-quiet
/tmp/ou3-export-quiet 3200 3200 > /tmp/ou3-quiet-word.jsonl
OPENBLAS_NUM_THREADS=1 python3 -m tools.stability.ou3_theorem.factor_diagnostic /tmp/ou3-quiet-word.jsonl --output /tmp/ou3-factor-diagnostic.json
python3 -m unittest discover -s tests/validation -p test_ou3_factor_word.py
```

The next experiment must vary complete admissible histories and carried
covariances while keeping factor/range structure. A positive diagnostic ratio
is only a feasibility gate. Uniform enclosure, full covariance upper bounds,
nonlinear remainder, finite-precision transfer and capture/retention remain.
