# OU-III stability proof handoff — PR #640

Branch: proof/ou3-imu-slow-fast

## Proof path

construction/capture -> magnetically informed H18 -> refinement/release -> regular A21 -> outer retention/finite inner entry -> local LaSalle strictness -> finite-error/prefix retention -> regime composition/float32.

Do not change estimator, tuner, calibration, runtime, gates, quality thresholds, MARINE bounds or MAGNETIC SERVICE to simplify proof.

## Physical assumptions now in force

IMU error is two-timescale for BOTH sensors:
e_a=b_a,s+b_a,f and e_g=b_g,s+b_g,f.
Slow components have amplitude/rate bounds and predecessor-linked increments. Fast components have amplitude plus temporal carried-window/primitive qualification; numerical H_a,C_a,H_g,C_g remain physically OPEN.

MARINE MOVING requires BOTH recurring gravity-direction span and displacement span:
Delta_g(W)>=theta_E on complete T_E windows and
Delta_p(W)=sup||p(s)-p(u)||>=P_E on complete T_P windows.
T_E,theta_E,T_P,P_E remain physically OPEN. Quiet water remains STILL. The old p=v=a=0 sin^3 rocking witness is historical only and is not MOVING.

Hypothetical proof-feasibility profile used only for mathematics:
T_E=T_P=H_a=H_g=60 s, theta_E=1 deg, P_E=.02 m,
C_a=1.2 m/s, C_g=.0042 rad.
It leaves physical ambiguity reserve 0.003056118 rad (~.1751 deg). Do NOT call these qualified device/vessel constants.

## Closed/reduced results on PR #640

* Local field-axis exclusion: inside sqrt(V)<=.15, physical force/field separation minus AW tube leaves 0.10527117647 m/s^2.
* Displacement excitation removes zero-translation MOVING aliases.
* Exact displacement boundary completion lambda(t)=T-t; q^T v is a LIN root coordinate.
* Held-H18 LIN structural detectability: separated S rows kill neutral v,p,S and OU AW root.
* Existing nuisance covariance upper proof lifts to held H18 first four LIN blocks.
* Held-H18 covariance/coefficient family compact; uniform homogeneous rho_L<1 exists qualitatively.
* Fixed-word affine source is bounded; all-time held-H18 LIN BIBO closes qualitatively.
* q^T v displacement boundary action and principal LIN release mean compactness close.
* Captured <=7 deg magnetic refinement gates have positive margin; MAGNETIC SERVICE gives finite captured-domain refinement/release. Compact A21 release set follows conditionally on capture.
* Zero homogeneous AG action does NOT imply zero base innovations; never revive that inference.
* Candidate 60-s profile passes physical ambiguity gate with 0.003056118 rad reserve.
* Old fast-gyro amplitude-only sinusoidal obstruction is excluded by candidate C_g temporal primitive.
* Causal terminal-AW reader action <=16 (norm <=4), but multiplying by outer radius is invalid.
* New paired_two_abel.py certificate layer pairs process/S and accelerometer physical-source columns BEFORE norms and telescopes internal slab endpoints exactly. It is NOT C_port,W.

## Controlling open bridge

Need

    sup_outer ||G_phys,W||_ind + E_nonlinear < 0.003056118 rad

for the hypothetical profile.

G_phys,W is the COMPLETE paired two-Abel physical-primitive -> normalized causal-source operator built from literal F,Q,K,H,S, AW sync, resets and the same physical p/v/S/a history.

Current estimator-word exporter has estimator factors/causal adjoints but does NOT export physical p/v/S/a lift columns aligned to every prediction/S/accelerometer boundary. Therefore a fixed-word G_phys,W cannot yet be formed from committed evidence.

Next executable task:
1. Extend the READ-ONLY stability observer/exporter to emit aligned physical p,v,S,a lift columns. Do not touch estimator/runtime behavior.
2. Construct G_phys,W using paired_two_abel.py, preserving process/accelerometer signs and collapsing internal boundaries before norms.
3. Run carried quiet/moving diagnostics only as feasibility evidence.
4. If feasible, build a source-uniform PSD/factor enclosure over the compact outer coefficient class. Do not promote replay extrema.
5. Add E_nonlinear/reset/arithmetic enclosure and compare with .003056118 rad.
6. If strict, return to literal six-column certificate: physical-to-nominal bridge => J_AG>0 => corrected-word full covariance upper and rho0<1 => outer entry.
7. Then evaluate linked chi_*/delta_ann and prefix retention.

## Other open obligations

* General construction/capture into the retained <=7 deg H18 domain.
* Magnetic-reference cone / Mahony proxy / continuous hard-iron transfer.
* Physical qualification of H_a,C_a,H_g,C_g,T_E,theta_E,T_P,P_E.
* Source-uniform corrected-word numerical loss/rho after J_AG.
* Outer A21 finite entry, linked finite-error/prefix supplies, regime composition and float32 closure.

## Forbidden/dead-end shortcuts

* Do not strengthen MARINE with a new EXCITED_MOVING assumption.
* Do not treat FAST as amplitude-only arbitrary residual.
* Do not set base innovations to zero from homogeneous zero action.
* Do not use pointwise global AW tracking.
* Do not use C_port,W as G_phys,W.
* Do not take separate process and accelerometer norms before pairing.
* Do not infer acceleration floor from P_E.
* Do not reverse covariance lower bounds into upper/action bounds.
* Do not promote carried rho, AW means, injection extrema or replay maxima.
* Do not return to obsolete O1/O2 kernel-ceiling architectures.
* Preserve coupled tau,sigma_aw,R_S,T_S chronology and actual event ordering.
