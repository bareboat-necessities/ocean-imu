# Coercivity of the actual finite binary32 covariance

## Result

For every **finite, exactly symmetric, strictly positive definite binary32**
21-by-21 covariance P, the prescribed joint24 storage matrix

`M(P) = diag(P^-1, I3)`

satisfies the universal mathematical bounds

`2^-133 I24 <= M(P) <= 2^5789 I24`.

The inverse and storage here are mathematical real objects, not a request to
compute these extreme constants or the inverse in binary32. No covariance
floor, eigenvalue clipping, source restriction or filter change is introduced.
The same bounds hold for the minimum-energy quotient metric and orthonormal
restrictions to unsupplied coordinates.

This proves the qualitative coercivity implication needed by the optional
history-wise tail theorem. A separate useful covariance envelope is no longer
required merely to establish that finite coercivity constants exist, **provided
all carried covariances stay finite, symmetric and SPD**. Establishing that
infinite-history membership is still an arithmetic/retention obligation.

## Exact proof

Every finite binary32 value is an integer multiple of q=2^-149 and has absolute
value at most F=(2^24-1)2^104. Hence P=q A for an integer matrix A. Strict
positive definiteness gives det A a positive integer, so

`det P = q^21 det A >= q^21`.

All eigenvalues are positive. Their sum is trace(P)<=21F=:B, so every
eigenvalue is <=B. Consequently

`lambda_min(P) >= q^21 / B^20`, and `lambda_max(P) <= B`.

Since B<2^133, we may weaken these to

`2^-5789 I21 <= P <= 2^133 I21`,

where 5789=149*21+133*20. Inverting the positive eigenvalues and appending I3
proves the claimed bound on M. This argument applies to the full carried
21-state covariance in both H18 and A21; held-bias cross terms are retained.
It is independent of source history, magnetic schedule and rounding mode,
as long as the resulting stored entries satisfy the stated format/domain.
It does not identify the actual recurrence with exact Riccati arithmetic.

For orthonormal complementary frames Q,N, the induced quotient energy is

`y^T M_q y = inf_z (Qy+Nz)^T M (Qy+Nz)`.

The lower bound follows from ||Qy+Nz||²=||y||²+||z||². For the upper bound,
choose z=0. Thus the same lower/upper constants hold for M_q. This is the
Schur-complement metric used by the regime diagnostic; it is not permission
to discard a non-invariant direction or to identify axial bias as a gauge.
An orthonormal coordinate injection Z likewise preserves the bounds on Z^T M Z.

## What this resolves and what it does not

The previous standalone demand for a history-wise numerical covariance upper
and lower envelope can be replaced, for this qualitative machine-storage
argument, by the existing finite-SPD domain preservation requirement. This
removes an independent coercivity-existence obstacle. A physically useful
ultimate radius still needs vastly sharper constants; the format bounds are
intentionally unsuitable as an engineering accuracy claim.

The following are **not** inferred:

- symmetry or positive diagonal implies SPD;
- real Joseph positivity implies rounded shipping positivity;
- finite covariance at observed endpoints implies indefinite finiteness;
- a below-one tangent ratio implies a nonlinear tail inequality;
- finite-format bounds establish source/target execution correspondence.

In particular, the shipping `symmetrize_Pext_()` averages mirrored entries;
it is not a positive-definiteness repair, and even its sum can overflow.
The compact covariance update uses subtractive arithmetic. Existing event
arithmetic, solve and rounded covariance obligations must establish finite SPD
membership at every required storage endpoint. Intermediate covariance/solve
sites also retain their existing totality requirements. Initial SPD alone
cannot be propagated by citing this format lemma.

The no-magnetic-service unipotent obstruction is unchanged. An ideal real
covariance can grow without bound; a finite binary32 implementation eventually
cannot follow that growth while preserving all other premises. A finite
subunit metric trace never proves an infinite strict contraction theorem.

## Exact native audit

`binary32_covariance_coercivity.py` checks exact binary32 representability,
exact symmetry and positive rational LDL pivots, then verifies determinant
separation on the dyadic grid. It does not accept a floating eigenvalue estimate
or repair a failed matrix.

The audit passes at all 36 retained native endpoints: root, sample 600 and
sample 1200 in each of twelve H18/A21/release and service/outage histories from
`service-regimes.json`. The report records the raw trace SHA-256, exact
positive determinant and minimum rational LDL pivot for each endpoint.
The slower-service histories remain outside-profile stress experiments.
The additional retained-service audit checks **18,586 recorded covariance
event boundaries** exactly: 6,195 in H18, 6,195 in A21 and 6,196 in the release
history, each including all 1,200 sample boundaries. Every matrix passes exact
finite-format, symmetry, rational positive-pivot and determinant-grid checks.
The report records operation counts and trace hashes. Internal partial writes,
pre-root covariance operations, moving sources and infinite prefixes are not
certified by these finite trace audits.

Reproduce after generating the service-regime native traces:

```sh
PYTHONPATH=.:tools/stability python3 -m \
  tools.stability.ou3_alt_contraction.binary32_covariance_coercivity \
  --native-directory /tmp/ou3-service-regimes \
  --audit-retained-traces --output /tmp/ou3-binary32-coercivity.json
```

The proof is a format-wide algebraic theorem. The native audit demonstrates
membership only at its recorded points; every ALT stability PASS remains false.
