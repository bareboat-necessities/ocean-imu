# H18 reference/AG/BG refinement and release compactness

There was a status/ledger inconsistency. The existing tail_stability.py contains a captured-domain proof of MagAutoTuner refinement liveness; it does not rely on a forced source timeout.

On captured tilt <=7 degrees, existing magnetic bounds |B| in [20,75] uT, horizontal component >=15 uT and measurement residual <=2 uT give positive margin at the literal 35% norm-variation and 5% horizontal-fraction gates. Hence every informative service sample is acceptable to the unweighted refinement tuner.

MAGNETIC SERVICE supplies at least one actually applied informative event per 1-s service window. The 128-sample refinement count therefore completes in finite uniformly bounded time; the 30-s elapsed-window gate is also finite. Refinement starts no earlier than 90 s. Conservatively charging all 250 internal accepted updates again plus the 1-s guard gives a finite captured-domain A21 release bound. This is a proof bound, not a runtime timeout.

AG/BG mean compactness during the finite bridge does not require the open general moving AG covariance theorem. Attitude lies on SO(3) and the retained capture chart; physical gyro bias is bounded and the implemented gyro-bias estimate is projected to .5 rad/s. LIN means are all-time BIBO. Hence all release means are bounded.

On this uniformly finite captured-domain interval, covariance, reference, frontend and tuner state evolve by finitely many continuous source operations on bounded inputs/coefficients. Their image is compact. Thus the A21 RELEASE SET is compact conditional on capture into the <=7-degree domain.

This does not prove general construction/capture into that domain and does not put release directly into sqrt(V)<=.15. The carried release audit already refutes that shortcut. The next theorem is outer A21 retention and finite inner entry from this compact release set, using temporal SLOW+FAST reachability, displacement/attitude excitation and actual MAGNETIC SERVICE.
