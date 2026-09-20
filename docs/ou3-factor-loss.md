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
