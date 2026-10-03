# Exact MOVING tilt--BA gauge obstruction

Status: **analytical real-arithmetic obstruction** to the proposed theorem that current qualified MOVING breaks the physical tilt/accelerometer-bias ambiguity. Failure class **B**. This is not a shipping-instability counterexample.

Use T_E=T_P=30 s, theta_E=2 deg, P_E=.03 m and set phi(t)=1 deg sin(2*pi*t/30), p(t)=.015 sin(2*pi*t/30)e_x, B=75e_x, delta=2 atan(1/200). Compare R0=Rx(phi) and R1=Rx(phi+delta), with zero gyro/FAST biases and b_a1=R0^T(a-ge_z)-R1^T(a-ge_z). Because a is parallel to e_x, the translational term cancels. The delivered accelerometer records coincide; the body rates coincide; B is invariant under both rotations, so magnetic records coincide. A common estimator initial state therefore produces the same complete deterministic shipping and accepted magnetic chronology.

The witness has max |p|=.015 m, max |v|=.003141592654 m/s, max |a|=.000657973627 m/s^2, max |jerk|=.000137805674 m/s^3, max |omega|=.003655409037 rad/s, |b_a1|=.098065274192 m/s^2, and max |dot b_a1|=.000358468690 m/s^3. Thus the SLOW BA amplitude/rate bounds hold and all FAST components are zero.

Every complete 30-s window contains a full period, so gravity-direction diameter is exactly 2 deg and displacement diameter is exactly .03 m. Hence both current MOVING clauses hold.

Therefore current MOVING + SLOW+FAST + MAGNETIC SERVICE cannot imply strict absolute physical tilt/BA kernel separation. Tightening interval arithmetic cannot repair that target. The viable theorem must use a measurement-compatible quotient/tube in BOTH STILL and MOVING, contract observable coordinates, and bound the gauge component without strengthening physical assumptions.

Executable arithmetic: tools/stability/ou3_theorem/moving_gauge_obstruction.py.

Structures preserved: actual 21-state shipping execution; current MARINE MOVING spans; one persistent SLOW+FAST split; magnetic-service chronology; Mahony, guard, tuner, scheduler, covariance/gains/Joseph, S=0 and BA projection.

Relaxations introduced: none.
