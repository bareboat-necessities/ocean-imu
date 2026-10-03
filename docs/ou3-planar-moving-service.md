# Planar MOVING symmetry reduction and magnetic-service target

Status: mean-coordinate symmetry proved analytically for the canonical real-arithmetic planar branch; all-time magnetic-service floor remains OPEN.

For the exact compatibility record use

    psi(t)=0.02 sin(pi t/10),
    y_g=psi_dot e_y,
    y_a=Ry(-psi)(a_x e_x-g cos(theta)e_z),
    y_m=75 Ry(-psi)e_x.

After exact symmetric magnetic acquisition/refinement the canonical reference is
B_ref=(75,0,0). On the planar branch the nominal world-to-body attitude is a
Y-axis rotation. The delivered accelerometer and magnetometer have zero Y
component. If AW, BA and LIN means lie in the XZ plane and gyro/BG lie on Y,
then prediction preserves that parity. Accelerometer and magnetometer residuals
remain in XZ. Their attitude Jacobians are -[v_xz]x, which map a Y-axis attitude
increment into XZ and do not create X/Z attitude increments from an XZ residual
on the parity-preserving covariance subspace. S=0 corrections use the same
block parity. BA radial projection preserves a zero Y component when BA is in
XZ. Thus the ideal real-arithmetic mean execution has a planar invariant
manifold. This explains the carried observations q_x=q_z=0 and hat b_a,y=0
without using them as proof premises.

The startup reference path is compatible with the same symmetry: yaw-stripped
tilt maps the planar magnetic record into the XZ plane, and MagAutoTuner
canonicalizes its reference to (horizontal magnitude,0,vertical component).
For the exact symmetric zero-hard-iron record the target reference is 75 e_x.
The continuous hard-iron estimator is exogenous and symmetry-compatible, but
its all-time eligibility/slew and the native float32 perturbation are separate
transfer obligations.

## Why the covariance does not collapse to one scalar

The mean symmetry does not make P diagonal. Pitch propagation uses the exact
attitude/BG transition, so the two heading/axial-BG probe pairs rotate/couple
inside the attitude/BG block. The magnetic innovation covariance is

    S_m = R_m + H_m P_theta,theta H_m^T,

and therefore depends on the carried covariance. A valid service proof cannot
replace S_m by a fixed noise matrix or independently choose gains.

The all-time certificate should therefore propagate only the reduced periodic
objects needed by magnetic service:

1. the attitude/BG covariance block and cross-blocks that enter its prediction;
2. any cross-covariances needed to reproduce P_theta,theta after literal
   accelerometer and S corrections;
3. the four homogeneous heading/axial-BG probe columns;
4. actual accepted magnetic events and their S_m;
5. the small set of tuner/scheduler variables that change those literal factors.

Everything is driven by the one 20-s periodic delivered record and the actual
startup state. This is a reduced HistoryCell, not a surrogate filter.

## Required inequality

For every placed one-second service window after the service root, certify

    lambda_min(sum Phi_i^T H_m,i^T S_m,i^-1 H_m,i Phi_i | E_hb) >= 1.

The finite carried diagnostic gives about 7.0248 on 102 disjoint windows, so the
available factor is about seven. That is feasibility evidence only. The rigorous
next step is a periodic/reachable enclosure of the reduced recurrence, preferably
a one-period invariant set plus a phase-uniform one-second information bound.

If the phase-uniform lower bound exceeds one and the reduced invariant set is
forward invariant, induction proves all-time MAGNETIC SERVICE on this history.
Then the already-proved identical-record pair and full-metric BA comparison imply

    max(V_+,V_-) >= 15.386492141376374 > 0.0225

at every common regular finite-SPD prefix, so universal retained absolute A21
entry is impossible on the fully admitted MOVING source class. That conclusion
must not be promoted before the service/invariance certificate closes.

Structures preserved: literal shipping mean/covariance/gain/Joseph/reset
chronology, startup/refinement/hard-iron path, S scheduler, coupled tuner,
MARINE MOTION, SLOW+FAST IMU and MAGNETIC SERVICE.

Relaxations introduced: the symmetry reduction removes coordinates only after
proving their parity invariance. No physical assumption or estimator behavior is
changed. Float32/all-time transfer remains open.
