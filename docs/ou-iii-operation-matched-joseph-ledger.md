# OU-III operation-matched Joseph/directional P4 route

This note records the proof route used after the first-accelerometer
sector-invariance diagnostic showed that shrinking the candidate angle cannot
close the deployed operation on the existing domain.

It is deliberately **not** a completed P4/P5 theorem claim.  It fixes the
acceptance calculus that the complete 18/21-state source-word backend must use.

## What is not changed

The stronger route preserves all deployment-facing assumptions:

- the P5 gauged entrance remains a **45 deg SO(3) geodesic** bound;
- the broad operation-matched geometry remains **0.80 rad / 45.84 deg**;
- the declared startup latent-acceleration error remains **0.3 g**;
- no `||delta a_w|| <= c sigma_applied` coupling is added to the operating
  domain;
- no filter gain, covariance equation, measurement schedule, state dimension,
  adaptation law, or runtime branch is changed.

The `30 -> 25 -> 20 -> 15` degree ladder remains useful as a numerical search
partition, but it is not a theorem assumption and it is not the closing route.

## Why per-operation sector invariance is the wrong gate

`ou3_p4_first_accel_sector_budget.py` asks whether the first accepted
accelerometer correction, considered in isolation, remains inside the 0.80-rad
outer attitude sector for every admitted source child.  On the current domain
its nuisance-over-budget ratio is above one at every candidate and remains
above one even at the 1-degree limit probe.

That result is a distance diagnostic, not a proof of instability.  The MEKF
uses a Joseph covariance update and an immediate left-error reset.  A transient
attitude-coordinate excursion is admissible if the complete source-word
Lyapunov/information level decreases.

Therefore the new promotion rule is:

> **Do not require each accepted operation to preserve the outer attitude
> sector.  Pair each operation with its own exact information change, preserve
> directional state blocks and source correlations, and test the complete
> return map.**

## Accelerometer: keep the effective `a_w` state coordinate

For the configured live accelerometer,

`J_aw = R_wb`

is orthogonal and full row rank.  If

`y_a = H_a z + eta_a`,

then with

`e_eta = R_wb^T eta_a`

one has exactly

`H_a E_aw e_eta = eta_a`

and hence

`K_a(H_a z + eta_a) = K_a H_a(z + E_aw e_eta)`.

The finite-angle residual is therefore not an unrelated measurement-noise term.
It is an effective source-correlated `a_w` state input.  In particular, the
large declared 0.3-g latent-acceleration error must **not** be paid as an
independent `eta^T R^-1 eta` penalty.

The latent-acceleration rotation itself is norm preserving.  What remains is to
propagate the exact effective state map through the source-correlated full
metric rather than multiply an independently maximized gain by an independently
maximized `a_w` error.

## Accepted measurements: use the Joseph information identity

For every accepted measurement,

`y = H z + eta`,

with Joseph posterior `P+`,

`z^T P^-1 z - (z-Ky)^T (P+)^-1 (z-Ky)`

`= y^T S^-1 y - eta^T R^-1 eta`.

The operation ledger specializes it as follows:

| operation | nonlinear treatment | information treatment |
|---|---|---|
| `S=0` | `eta=0` exactly | direct positive Joseph decrease |
| accelerometer | finite-angle residual absorbed into effective `a_w` state coordinate | Joseph decrease on the joint effective state; no standalone vector-eta penalty |
| magnetometer | radial residual has exactly zero Kalman-gain action; useful residual is an effective tangent coordinate | Joseph decrease on the effective tangent coordinate |
| rejected/not due | exact identity | zero change |
| quaternion injection/reset | exact covariance congruence | no condition-number multiplier; retain explicit Cayley reset defect `rho` |

The finite-angle accepted accelerometer/magnetometer packet already has a strict
positive directional information lower bound from
`ou3_p5_outer_information_geometry.py`.  The remaining task is not to discover
positivity; it is to prove that prediction, effective-state transport, reset
defects, and full-state cross terms do not consume that decrease on any
source-reachable complete word.

## CI gate added by this route

`ou3_p4_operation_matched_joseph_ledger.py` fails if any of the following is
introduced:

1. the declared 45-degree P5 entrance is shrunk;
2. the `a_w/sigma_applied` consistency condition is silently promoted into the
   domain;
3. the old per-operation sector-invariance test is restored as a P4 promotion
   gate;
4. the accelerometer finite-angle residual is reclassified as independent
   measurement eta;
5. the magnetometer radial residual is charged despite exact gain annihilation;
6. a reset covariance condition-number multiplier is reintroduced; or
7. this intermediate ledger is relabeled as complete P4/P5 closure.

## Remaining numerical closure

The next backend must propagate, on the same recurrent source paths and with
outward rounding:

1. the source-correlated effective-state map `z -> z_eff`;
2. each accepted operation's Joseph information decrease;
3. the exact Cayley/quaternion reset defect `rho`;
4. exact finite-angle prediction defects;
5. the complete 18-state H and 21-state A covariance/information cross terms;
6. the normalized translation/nontranslation cross block below the existing
   lower-enclosed Schur budget; and
7. the complete return-map decrease while allowing transient attitude-sector
   excursions.

Only after that complete-word decrease is established may P4 be promoted and P5
compute a finite capture count to the inner stochastic localization level.
