# Bias-gauge boundary from the two-timescale IMU model

For an exact zero-translation alias producing the same level accelerometer packet with zero FAST accelerometer error, let d=R^T e_z. The required SLOW bias is b_a,s=g(d-e_z)+c, where the common representation offset c cancels in differences. If a complete MARINE window contains gravity directions separated by theta, then

    |b_a,s(t2)-b_a,s(t1)| = 2 g sin(theta/2)
                           <= min(2 B_a,s, D_a,s T_E).

Zero measured gyro with zero FAST gyro error also gives omega=-b_g,s. Gravity-direction path length is no larger than integral |omega|, hence theta <= B_g,s T_E. Therefore every waveform in this slow-only quiet-alias class must satisfy

    theta_E <= Theta_slow(T_E)
    Theta_slow = min(2 asin(min(2 B_a,s,D_a,s T_E)/(2g)),
                     B_g,s T_E, pi).

The inherited candidate budgets give Theta_slow=0.00173352 rad (0.09932 deg) at 17 s, 0.00611831 rad (0.35055 deg) at 60 s, and 0.01019721 rad (0.58426 deg) at 100 s. Thus a QUALIFIED existing MARINE pair T_E=60 s, theta_E>0.35055 deg would exclude every all-SLOW quiet alias of this class. This does not adopt that pair: current T_E and theta_E remain OPEN.

The existing phi=.001 sin^3(t/40) construction remains below this boundary for its symbolic T_E=80 pi, theta_E=.002, so a theorem over arbitrary positive symbolic excitation cannot be closed from SLOW bounds alone.

For two complete delivered endpoint cells, a qualified FAST profile extends the necessary condition to

    2 g sin(theta_E/2)
      <= min(2 B_a,s,D_a,s T_E) + c_a(d0)+c_a(d1),

where c_a(d)=min(B_a,f,C_a/min(H_a,d)). The strict reverse inequality excludes this alias. If H_a,C_a are OPEN, this conclusion is OPEN. Using 2 B_a,f as an admission rule would recreate the old loophole.

This endpoint theorem is not promoted to full gauge exclusion. A general pair of physical histories also requires two reachable error histories, translation, transported gyro compatibility and actual MAGNETIC SERVICE. Those constraints can only shrink the set, but their source-uniform intersection remains to be proved.

Consequences: the old arbitrary-residual witness is no longer automatically admitted; .01 sin(.5t) cannot be wholly SLOW and requires a large FAST accumulation; the smaller sin-cubed witness is genuinely SLOW and survives weak MARINE pairs. No stronger MARINE assumption is introduced. The existing T_E,theta_E and both FAST profiles must be qualified before the remaining joint compatibility obstruction can be decided.
