# Held-H18 LIN compactness and BIBO

The displacement boundary row q^T v is a coordinate of z=(v,p,S,a_w); no second reader is introduced.

Three distinct actual S observations kill the neutral (v,p,S) root because rows [t^2/2,t,1] have determinant (t1-t0)(t2-t0)(t2-t1)/2 != 0. A fourth distinct S observation kills the integrated-OU a_w extension by the existing extended-Chebyshev/four-S lemma. Prediction and covariance-matched corrections are nonexpansive in the literal covariance metric.

The covariance upper bound was already proved in docs/ou3-nuisance-upper-proof.md. After 17 s it bounds P_nn using a dominating comparison that keeps autonomous LIN prediction and actual S corrections, fixes realized AW sync increments, and OMITS accelerometer and magnetic corrections. Its cancellation of an arbitrarily large neutral root uses three separated S observations. BA is a separate fifth-block argument.

Therefore the first four blocks do not require active BA. The same tuner, pseudo-S scheduler and LIN prediction operate during regular Live H18 before BA release. Holding BA does not alter LIN F/Q or suppress S corrections. Hence

    P_LL <= 5 diag(b_v^2 I,b_p^2 I,b_S^2 I,156^2 I)

at every operation boundary after 17 s of regular held-H18 service, with LIN cross covariance retained. Other optimal corrections only decrease the LIN principal covariance.

The same source proof supplies dt in [.004,.006], tau in [.02,12], sigma_aw<=4, R_S sigma<=100 and applied-S gap <=.156 s. Positive measurement-noise floors make gains continuous on this bounded family. Thus held-H18 homogeneous words after 17 s form a compact coefficient/covariance family.

Compactness plus the zero-action kernel gives by continuity a source-uniform strict homogeneous factor on a fixed separated four-S word:

    V_L,end <= rho_L V_L,root, rho_L<1.

No numerical rho_L is claimed. The remaining BIBO step is AFFINE forcing: AG/BG and held BA enter through literal cross-coupled gains; physical SLOW+FAST sensor histories and physical S enter base innovations. A fixed word is finite, but this source vector must be bounded on the SAME physical history.

Once a source-uniform fixed-word input bound d_H exists,

    sqrt(V_L,end) <= sqrt(rho_L) sqrt(V_L,root)+c_H ||d_H||,

and iteration gives an all-time H18 LIN radius. The covariance upper bound then converts storage to a finite Euclidean velocity bound, closing q^T v, the displacement boundary row and principal LIN release compactness.

Thus Riccati compactness is no longer the blocker. The next obstruction is the linked affine H18 input bound under MARINE + SLOW/FAST IMU + MAGNETIC SERVICE.
