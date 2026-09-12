# OU-III parallel ALT proof handover

## Resume point and independent scope

Continue PR #523 on its current branch, preserving concurrent commits. Read
`AGENTS.md`, `docs/ou3-alt-proof-plan.md`, this handover, the ALT research ledger,
`docs/ou3-alt-finite-measurement-proof.md`,
`docs/ou3-alt-finite-core-composition.md`, `docs/ou3-alt-runtime-primitives.md`,
`docs/ou3-alt-mahony-binary32.md`, and `docs/ou3-brmm-main-handover.md`.
The original P2/P3/P4/P5 route remains independently continuable; its premises
and gates must not be weakened. `P3=1e-18` remains frozen.

## Certified ALT deployment scope

ALT explicitly excludes the optional wind-heel retarget feature. Certified
histories require `wind_heel_rad_ == 0` from construction onward and zero calls
to `update_wind_heel()`. Shipping is unchanged and already defaults heel to
zero. Hence B'=B throughout the certified history, de-heeling is identity, and
no wind-heel/body-frame retarget event belongs to the ALT hybrid language.
`deployment_scope.py` and the storage guards fail closed if this scope is lost.

Retain joint24 `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)`, full 21-state covariance,
all motion/bias cross terms, one persistent Live S origin and corrected
COMPLETE-BRMM including acceleration <=8.8 m/s^2 and all-time centered primitive
`D_S<=1100 m*s`. Preserve same-signal WPE/bandpass/sigma/tau/T_S ancestry,
staged tuner state and separate BIAS0/1/2 histories.

## Finite runtime graph already available

Local finite identities retain actual continuous physical attitude increments
and angular/model defect, correlated translation moments, one shared physical
accelerometer-bias driver, inverse-free `K Sigma=N`, accepted/rejected SafeLDLT
branches, quaternion injection, same-beta projection and full Joseph/reset
covariance. H18 retains latent BA covariance.

Prediction/runtime modules now generate the ordinary shipping attitude,
OU-chain, Qaxis and BA covariance coefficients from same-history runtime roots;
pending a_w floor, S scheduler/service, accelerometer vibration guard, applied
R_acc, held-sample accelerometer forcing, staged tuner commits and the first
ordinary Live IMU prefix are represented. The initialized private Mahony path
has a named exact binary32 graph under its stated arithmetic profile. Target
compiler/libm/profile qualification remains open.

## Startup magnetic source and capture

The theorem-facing startup magnetic path is structurally attached through the
persistent Mahony proxy, world-frame gravity gate, asynchronous admission clock,
default MagAutoTuner recurrence, yaw-stripped tilt frame, physical magnetic
source and the same `PhysicalKinematics` endpoint as the main word.

`MAG-BMM150-DET-v1` admits commissioned installations with

- `20 <= ||B_W|| <= 75 uT`;
- horizontal field `>=15 uT`;
- body hard iron `<=5 uT`;
- deterministic residual `<=2 uT` per theorem sample.

Together with the declared startup gravity-direction error <=0.02 rad, the
real-arithmetic startup proof gives total horizontal perturbation <=8.5 uT,
`|sin(delta_yaw)| <= 17/30`, yaw error <0.61 rad and total attitude error
<0.63 rad <45 deg. No statistical `1/sqrt(N)` reduction is used.

The generation-zero magnetic reference write and pending absolute-yaw gauge are
kept separate from the later MEKF handoff. The zero-heel handoff installs the
shipping tilt/yaw covariance and retains the accepted handoff quaternion in the
MEKF's internal W->B convention.

## Fresh H18 entry no longer uses an assumed error box

`finite_startup_live_entry.py` composes shipping `goLive`/`enterLive_`: it seats
P_aw,aw on the same committed Sigma_aw, clears every a_w cross-covariance,
requires committed Live R_S and retains accelerometer-bias learning disabled.

`finite_fresh_joint24_entry.py` then derives the actual fresh joint24 state from
the SAME physical `Reference` and estimator coordinates:

`c = Cayley(q_true_WB * conjugate(q_hat_WB))`,
`e_bg=b_g- b_g_hat`, `e_v=v-v_hat`, `e_p=p-p_hat`,
`e_S=S_centered-S_hat`, `e_aw=a-a_w_hat`,
`e_ba=beta-b_a_hat`, and the final three coordinates are the same `beta`.

The one-time Live origin and fresh centered physical S=0 are enforced. The
result is the existing full-covariance H-mode `finite_core.State`. Covariance
consistency and an independently chosen fresh-entry radius are not premises.
Deployment binary32 correspondence and later basin membership are still open.

## H18/A21 hybrid control

`MAG-CALL-SCHEDULE-v1` requires the first post-Live magnetometer call within
40 ms and subsequent gaps <=40 ms. Because shipping counts attempted post-delay
`updateMag()` calls independently of innovation acceptance, the 250-count unlock
guard clears within 10 s of gauged Live and the strict >1 s condition is then
automatically satisfied.

Do NOT assume eventual A21 under an arbitrary external accelerometer-bias hold.
`finite_h18_a21_transition.py` retains both literal continuations:

- no hold: the qualified unlock edge applies the exact H18->A21 BA variance
  floor;
- hold active: the internal lock clears but the filter may remain H18
  indefinitely;
- first later hold release while Live applies the same exact H18->A21 edge;
- asserting a hold from A21 returns to H18 and zeros BA cross-covariances exactly
  as shipping does.

Thus both H18 and A21 continuations belong to the theorem's hybrid language.

## Decisive remaining object

**The complete source-uniform finite 600-step runtime word is still open.**
The current limiter is no longer wind heel or a missing fresh-entry error set.
The main remaining work is:

1. deployment/binary32 correspondence for startup `atan2`, AngleAxis,
   quaternion normalization, handoff setters and clocks;
2. connect the exact fresh H18 `CORE.State` to the first represented Live IMU,
   pseudo-S and asynchronous-mag prefixes with no frontend/tuner/scheduler gap;
3. attach the complete qualified asynchronous magnetometer/refinement schedule
   and remaining continuous hard-iron branches required by the declared scope;
4. close remaining solver/hygiene/finite-precision branches and source-uniform
   COMPLETE-BRMM/BIAS ancestry across every literal prefix;
5. only then allow `assert_finite_storage_master` to admit a high-precision
   feasibility diagnostic/common joint24 storage search.

Do not substitute traces, random seeds, frozen gains, independent coefficient
boxes, covariance consistency or a convenient fresh-entry set for these tasks.

## Gate state

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
Original P4/P5 remain independent and unpromoted. No certified rho, retained
basin or ultimate bound is claimed yet.
