# OU-III proof research state

This file is the current proof ledger. Historical routes belong in Git history.

## Active theorem

The target is **uniform regional practical stability with finite capture for
perturbations of physically admissible multimodal directional JONSWAP SEA3 sea
states**.  P3 is a recurrent finite-window linear certificate; P4 is nonlinear
finite-window dissipation; P5 is finite capture.

## Canonical P3 architecture

P3 strictness is now assigned to the strongest corrective mechanism actually
visible in the shipping filter: recurrent measurement innovation dissipation,
with the S=0 pseudo-measurement and its adaptive R_S schedule as the primary
translation correction.

For every exact linear Kalman measurement update,

`V^- - V^+ = r^T S_innov^-1 r`, `V=e^T P^-1 e`.

Prediction is non-expansive because process covariance is PSD.  Summed over a
recurrent word, the sequential innovations are the block-Cholesky whitening of
the batch measurement record, so

`V_0 - V_W >= e_0^T D_W e_0`,

`D_W = O_W^T Sigma_Y^-1 O_W`.

The source-uniform process UCC certificate supplies a covariance lower `L_W <=
P_0`; it is a metric scale, not the primary strictness mechanism.  The P3 gate
is the full-matrix inequality

`D_W >= delta L_W^-1`, `delta >= 1e-18`,

which implies `D_W >= delta P_0^-1` and hence

`V_W <= (1-delta)V_0`.

### Why R_S is central

The shipping S=0 update has residual `-delta_S`, innovation covariance
`P_SS+R_S`, and gain built from the complete cross-covariance `P(:,S)`.  It
therefore corrects the correlated `v,p,S,a_w` chain rather than only the S
coordinate.  Three separated recurrent S observations reconstruct the endpoint
`v,p,S` integrator coordinates; the stable a_w direction and the required
Normal-Live accelerometer updates close the fourth translation direction.

R_S and pseudo cadence may not be replaced by independent extrema.  The
shipping tuner couples tau, sigma, R_S and T_S at one operating point.  The
Cubic law uses the tau^3 drift-band schedule with cadence information-rate
normalization; the Riccati and SpectralMSE laws explicitly include the realized
T_S.  The proof must preserve that same-source coupling.  A Cartesian corner
such as maximum R_S combined with unrelated minimum sea/process scale is
inadmissible unless the actual SEA3+tuner dynamics permit it.

### SEA3 preconditions consumed by P3

The recurrent word uses the declared Normal-Live conditions: every valid IMU
sample executes the accelerometer update, accelerometer rejection is outside the
theorem branch, vector PE recurs within the declared window, lever arm is off,
and the active vibration-guard branch is outside the current proof scope.
Physical SEA3 height/period and partition-energy coupling is retained rather
than replaced by independent rectangular extrema.

### Routes explicitly retired

The canonical P3 proof must not use:

- one-sample strict Riccati process injection;
- per-sample SPD covariance lowers;
- commit-aligned source-word propagation;
- selected-process-mode strictness as the primary contraction mechanism;
- determinant/trace eigenvalue scalarization;
- one scalar information-beta attenuation;
- the retired 800-state P2 graph or predecessor/source-history enumeration.

The endpoint covariance upper remains useful for boundedness and P4, but the
old one-step `~2.12e-35` number is diagnostic only and cannot gate P3.

## Current quantitative obligation

There are exactly two open P3 matrices:

1. `D_W`: the SEA3-coupled R_S/T_S weighted batch innovation-information lower,
   including recurrent S=0, required accelerometer packets, and vector PE;
2. `L_W`: the source-uniform UCC covariance lower in the same H18/A21
   coordinates.

After those are emitted, run one validated full-matrix comparison
`D_W >= delta L_W^-1` for H18 and A21 with the unchanged `1e-18` gate.  P4
remains blocked until that passes; P5 remains blocked until P4 closes.
