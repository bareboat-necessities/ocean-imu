# Literal six-column J: joint AW/injection obstruction

For the hypothetical profile T_E=T_P=H_a=H_g=60 s, theta_E=1 deg,
P_E=.02 m, C_a=1.2 m/s, C_g=.0042 rad, the physical ambiguity gate has
positive margin. The remaining G0 premises are estimator-side.

Theorem G0 already proves an injection-free aggregate six-column floor
s^2>=1.486786e-3 if two 16-s nominal statistics satisfy
m_perp<=.4 m/s^2 and mean|f_hat|/g<=1.2. Carried words satisfy .348 and 1.091,
but this is not source-uniform.

The literal AW loop has an exact signed identity. Physical acceleration
increments telescope; the nominal 16-s mean is a signed gain-weighted sum of
AW corrections plus root/OU terms. Absolute correction accumulation is a
documented kill test and cannot prove .4. The candidate FAST temporal caps
bound the sensor contribution to this signed sum, but attitude/held-BA terms
share the same correction gains and outer state.

The literal injection transport likewise has an exact signed group identity:
angle(J)<=angle(M0)+angle(MN)+integral|omega_tilde|. The old perturbative
sum of injection norms is a documented dead end. Candidate C_g controls the
FAST gyro integral, while slow gyro level/attitude evolution must be controlled
jointly by actual MAGNETIC SERVICE (Lemma T), not by the .5-rad/s estimator
invariant alone.

Therefore the two apparent scalar premises cannot be certified independently
by separate worst-case boxes. A source-uniform proof must use a JOINT
zero-loss/compactness contradiction on the compact outer class: assume a
sequence of literal corrected words with sigma_min(AG rows)->0. Compactness
gives a limiting same-history word. G0/Lemma T force its nominal AW signed
means toward the collinear compatibility class; the exact AW loop identity,
MARINE displacement+attitude excitation and temporal SLOW+FAST constraints
then restrict the physical forcing; the signed injection identity plus
MAGNETIC SERVICE restricts the rotating kernel. If the only limit is the
already-excluded zero-action compatibility trajectory, compactness yields
a positive literal six-column floor without separately proving .4 and an
injection norm threshold.

This is now the shortest viable route. The repository does not yet contain
that joint contradiction. Claiming the two scalar inequalities independently
from the candidate profile would be false. No estimator or physical assumption
needs to change; the remaining task is one compactness/nullspace lemma on the
literal six-column array.
