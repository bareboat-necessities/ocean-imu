# Held-H18 LIN BIBO theorem

The displacement boundary row q^T v is a coordinate of z=(v,p,S,a_w); no second reader is introduced.

Three distinct actual S observations kill the neutral (v,p,S) root, and a fourth kills the integrated-OU a_w extension by the existing extended-Chebyshev/four-S lemma. Prediction and covariance-matched corrections are nonexpansive in the literal covariance metric.

The existing nuisance-upper proof lifts to held H18 because its first four blocks use only autonomous LIN prediction, actual S corrections and fixed realized AW sync increments while omitting accelerometer/magnetic corrections. BA activity appears only in its separate fifth block. Thus after 17 s of regular held-H18 service,

    P_LL <= 5 diag(b_v^2 I,b_p^2 I,b_S^2 I,156^2 I)

at every operation boundary, with LIN cross covariance retained. The same proof supplies compact dt/tau/sigma/R_S/scheduler ranges. Positive innovation-noise floors make actual gains continuous. Compactness plus the zero-action kernel gives a source-uniform homogeneous rho_L<1 on a fixed separated four-S word.

For BIBO, unlike the later sharp outer-entry inequality, only FINITENESS of the affine word forcing is required. On that fixed word dt>=.004 gives finitely many operations. MARINE supplies finite p,v,a,jerk,angular-rate and primitive bounds. SLOW accel/gyro biases have finite amplitudes; FAST errors have finite amplitudes. The held BA error is bounded after a feasible completed projection, and the implemented gyro estimate is bounded. Physical S is bounded by the persistent primitive. Gravity/magnetic/reference envelopes are finite. Since covariance/tuner/scheduler and actual gains form a compact family, the finite chronological affine word map is continuous on a compact input/coefficient product. Therefore a source-uniform constant D_H<infinity exists with

    ||d_H||_word <= D_H.

This existence proof does NOT use independence, zero mean, or an unrestricted residual channel. The OPEN fast temporal H,C values are not needed for mere BIBO because the fixed word has finite event count and B_f is finite. They remain essential for the later same-history gauge/outer-entry margin, where replacing temporal reachability by independent amplitude boxes would be too weak.

Hence

    sqrt(V_L,k+N) <= sqrt(rho_L) sqrt(V_L,k) + c_H D_H

and iteration gives a finite all-time held-H18 LIN radius after the 17-s entry prefix. The finite prefix is bounded by continuity from the finite construction/handoff state and bounded inputs. The covariance upper comparison converts the storage radius to finite Euclidean v,p,S,a_w bounds.

Consequences now CLOSED at the qualitative/source-uniform existence level:
* all-time held-H18 LIN BIBO;
* q^T v displacement boundary action;
* principal LIN mean compactness at A21 release.

No numerical BIBO radius is claimed. This does not close H18 reference-refinement time, AG/BG compactness, nonlinear retention, or the sharp linked outer-to-inner A21 supply. For those later steps, retain the SLOW+FAST temporal reachable sets and actual MAGNETIC SERVICE rather than this deliberately coarse amplitude compactness argument.
