# OU-III operation-matched Joseph/directional P4 route

This note records the stronger P4 route after the first-accelerometer
sector-invariance diagnostic showed that shrinking the candidate angle cannot
close the deployed operation on the existing domain.

It is deliberately **not** a completed P4/P5 theorem claim. It fixes the
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

The `30 -> 25 -> 20 -> 15` degree ladder remains useful as a diagnostic/search
partition, but it is not a theorem assumption and it is not the closing route.

## Why per-operation sector invariance is the wrong gate

`ou3_p4_first_accel_sector_budget.py` asks whether the first accepted
accelerometer correction, considered in isolation, remains inside the 0.80-rad
outer attitude sector for every admitted source child. On the current domain
its nuisance-over-budget ratio is above one at every candidate and remains
above one even at the 1-degree limit probe.

That result is a distance diagnostic, not a proof of instability. The MEKF
uses a Joseph covariance update and an immediate left-error reset. A transient
attitude-coordinate excursion is admissible if the complete source-word
Lyapunov/information level decreases.

Therefore the new promotion rule is:

> **Do not require each accepted operation to preserve the outer attitude
> sector. Preserve each operation's directional information and signed
> nonlinear term, transport both through the exact return map, and test the
> complete source word before reducing to a scalar margin.**

## Exact structural rank: why the word must be directional

`ou3_p4_directional_packet_rank.py` gives a structural, not numerical, result.
In the rotation gauge used by the vector proof,

`H_a = [-[f]_x, I_aw]`

and

`H_m = [-[m]_x, 0_aw]`.

The accelerometer has exact rank three because the `a_w` block is full rank.
The magnetometer skew map has exact rank two for the admitted nonzero magnetic
vector. The stacked vector packet has the exact nonzero null family

`delta theta = alpha m`,

`delta a_w = [f]_x delta theta`,

so its rank is exactly **5** on the 6-dimensional H-mode `(theta,a_w)` block.
In A mode the same five measurement directions act on the 9-dimensional
`(theta,a_w,b_a)` block, giving nullity four.

Thus an instantaneous strictly positive full-state measurement-information
margin is algebraically impossible. Even with a due `S=0` pseudo measurement,
the direct same-sample rank is only eight. P3 achieves global detectability
because prediction transports these null directions and recurrent directional
information accumulates over the word. P4 must preserve the same mechanism.

## Accelerometer: exact correction-range reduction, not eta deletion

For the configured live accelerometer,

`J_aw = R_wb`

is orthogonal and full row rank. If

`y_a = H_a z + eta_a`,

set

`e_eta = R_wb^T eta_a`.

Then exactly

`H_a E_aw e_eta = eta_a`

and

`K_a(H_a z + eta_a) = K_a H_a(z + E_aw e_eta)`.

This is an exact reduction of the **state correction** to an effective
source-correlated `a_w` coordinate. It means the large declared 0.3-g
latent-acceleration state error is not an unrelated measurement disturbance and
must not be multiplied by an independently maximized gain as in the retired
sector-invariance diagnostic.

It does **not** mean that nonlinear measurement `eta_a` becomes zero in the
Lyapunov calculation.

## Accepted measurements: preserve the signed Joseph identity

For every accepted measurement,

`y = H z + eta`,

with Joseph posterior `P+`, the exact identity is

`z^T P^-1 z - (z-Ky)^T (P+)^-1 (z-Ky)`

`= y^T S^-1 y - eta^T R^-1 eta`.

The stronger route keeps both terms together on the same source cell. In
particular, it does **not**:

- drop `eta^T R^-1 eta` after introducing an effective coordinate;
- replace it by a separately maximized `||eta||` bound; or
- turn one directional packet into a fictitious positive scalar full-state
  margin.

The operation ledger specializes the identity as follows:

| operation | nonlinear treatment | information treatment |
|---|---|---|
| `S=0` | `eta=0` exactly | direct PSD directional credit `y^T S^-1 y` |
| accelerometer | finite-angle residual represented exactly in the effective `a_w` correction range | retain `y^T S^-1 y - eta_a^T R_a^-1 eta_a` jointly; never maximize eta independently |
| magnetometer | radial residual has exactly zero Kalman-gain action; useful residual is an effective tangent coordinate | cancel/retain the exact source-correlated tangent Joseph form; no independent radial penalty |
| rejected/not due | exact identity | zero change |
| quaternion injection/reset | exact covariance congruence | no condition-number multiplier; retain explicit Cayley reset defect `rho` |

`ou3_p5_outer_information_geometry.py` still supplies a strict positive
**attitude-geometry** vector-pair constant on the finite-angle handoff nodes.
That constant is useful geometry, but it cannot be promoted to a full-state
packet credit because the rank-five theorem proves that such an instantaneous
scalar credit cannot exist.

## CI gate added by this route

The new producers fail if any of the following is introduced:

1. the declared 45-degree P5 entrance is shrunk;
2. the `a_w/sigma_applied` consistency condition is silently promoted into the
   domain;
3. the old per-operation sector-invariance test is restored as a P4 promotion
   gate;
4. the 0.3-g `a_w` state error is reclassified as independent measurement eta;
5. the finite-angle Joseph eta term is dropped or independently maximized;
6. a rank-five packet is relabeled as a positive scalar full-state packet
   margin;
7. the magnetometer radial residual is charged despite exact gain annihilation;
8. a reset covariance condition-number multiplier is reintroduced; or
9. this intermediate calculus/rank result is relabeled as complete P4/P5
   closure.

## Remaining numerical closure

The next backend must, on the same recurrent source paths and with outward
rounding:

1. build the PSD directional form for each accepted `S`, accelerometer and
   magnetometer operation and the associated signed nonlinear eta form;
2. transport those forms through source-correlated prediction, effective-state
   maps and every exact quaternion/reset congruence;
3. accumulate the directional forms over the complete recurrent H/A word
   **before** taking any scalar lower eigenvalue;
4. enclose the exact Cayley/quaternion reset defect `rho` and finite-angle
   prediction defects on those same cells;
5. retain all H=18 / A=21 covariance/information cross terms;
6. outward-enclose the normalized translation/nontranslation cross block below
   the existing lower-enclosed Schur budget; and
7. prove a strict complete return-map decrease while allowing transient
   attitude-sector excursions.

Only after that complete-word decrease is established may P4 be promoted and P5
compute a finite capture count to the inner stochastic localization level.
