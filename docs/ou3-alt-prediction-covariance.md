# ALT finite prediction-covariance block lemma

## Scope

This lemma removes the previous arbitrary full `F21,Q21` description of the
ordinary OU-III covariance prediction. It covers the shipping state order
`[dtheta,b_g | v,p,S,a_w | b_a]` at the point after the ordinary prediction
blocks and before pending `a_w` covariance inflation, final symmetry hygiene,
and periodic S=0 service. Those later events remain separate literal prefixes.
It does not change the deployed filter or authorize a storage search.

## Literal block relation

For the current shipping configuration with gyro and accelerometer bias states,
write `F_AA,Q_AA` for the actual six-state attitude/gyro-bias prediction,
`F_LL,Q_LL` for the 12-state integrated-OU block, and `phi_b,Q_BB` for the
accelerometer-bias branch. Before the later post-prediction events, the code is
exactly the block-diagonal congruence

```
F = diag(F_AA, F_LL, phi_b I3),
Q = diag(Q_AA, Q_LL, Q_BB),
P+ = F P F' + Q.
```

This statement retains every AA/LL/BA cross covariance: `P_AL` receives
`F_AA P_AL F_LL'`, while both BA cross blocks use the SAME `phi_b`. H18 is the
literal held branch `phi_b=1,Q_BB=0`; A21 uses its active bias factor and process
covariance. No 18-state marginal covariance is introduced.

`finite_prediction_covariance.py` evaluates both the shipping block-update order
and the dense congruence in exact rational arithmetic and requires entrywise
equality. The equality is polynomial in the block entries and predecessor P, so
its proof is not a sampled stability experiment.

## Same mean coefficients

`F_LL` is no longer a free matrix. The constructor consumes the SAME per-axis
`(phi_va,phi_pa,phi_Sa,alpha,h)` tuple already used by the finite mean predictor:

```
[1    0 0 phi_va]
[h    1 0 phi_pa]
[h^2/2 h 1 phi_Sa]
[0    0 0 alpha]
```

Thus the mean and covariance linear transitions cannot silently use unrelated
OU coefficients. For correlated `a_w`, `Q_LL` is assembled exactly as shipping:
each 3x3 group block is `Sigma_aw * Qaxis_unit[g,h]`. The independent-axis
fallback is also represented without cross-axis terms.

## What remains conditional

This closes a structural covariance-composition gap, not the runtime/source
graph. The actual source relation still has to generate `F_AA,Q_AA`, the
analytic `Qaxis` values, active BA `phi_b/Q_BB`, and the same-history tuner
quantities from the runtime predecessor. The exact-attitude-Q branch includes
trigonometric/Simpson and PSD-hygiene decisions; these are not replaced by free
noise rectangles. Pending `a_w` inflation, symmetry/nonfinite repairs and
periodic S=0 scheduling remain separate events that must be composed literally.
Finite-precision effects remain open.

Consequently the source-uniform 600-step word, common storage, rho, retained
basin and startup theorem remain unproved. All ALT gates remain false. The
original P2/P3/P4/P5 proof route is untouched.
