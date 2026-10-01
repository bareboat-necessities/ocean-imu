# Literal A21 release outer-storage audit and source-uniform entry result

This note performs the two requested calculations after the radius-local
field-axis lemma. It does not promote a finite replay to a theorem.

## 1. Literal BA-eliminated release storage

The dedicated unchanged-header diagnostic snapshots the first operation at
which shipping reports acc_bias_updates_enabled(). On the carried diagonal-wave
construction this occurs at step 24016, t=120.079997316... s.

For storage V=e'P^-1 e, eliminating BA means

    V_o(e_o)=min_z [e_o;z]'P^-1[e_o;z]=e_o'P_oo^-1 e_o.       (RE1)

At the release boundary P_ob=0 because held BA was explicitly decoupled, but
RE1 remains the correct elimination identity after cross covariance regrows;
the covariance Schur complement is a conditional, not eliminated, covariance.

Using the literal release covariance/state and the exact physical construction

    p=-(3/2)sin(2t)(1,0,1),
    v=-3cos(2t)(1,0,1),
    a=6sin(2t)(1,0,1),

with R=I and zero true biases, the carried release gives

    V_o = 7.9300039948e4,
    sqrt(V_o) = 281.6026.                                    (RE2)

The target field-alignment ball is r_FA=.15, r_FA^2=.0225, so this finite
release is very far outside it.

To identify the obstruction without attributing covariance cross terms to
arbitrary coordinates, minimize the outer quadratic over every coordinate
except each named principal block. The resulting rigorous block lower bounds
on the same carried release are approximately

    block       V_block,min       sqrt(V_block,min)
    attitude       9.7976              3.1301
    gyro bias      3.59e-7             5.99e-4
    velocity     331.982              18.220
    position    1996.290              44.680
    S             84.977               9.218
    AW          60867.918             246.714.               (RE3)

Thus the dominant non-BA obstruction is AW, not attitude. The literal physical
AW error at this release is approximately

    e_aw=(-5.9116273,-0.0001332,-5.9459625) m/s^2,

norm about 8.385 m/s^2. Attitude alone also excludes .15 on this replay:
the physical tilt is only .39342 deg, but the attitude covariance eigenvalues
are already 3.046e-6--6.154e-6, giving V_theta,min=9.7976.

RE2--RE3 are finite carried diagnostics only. The construction is not an
all-time MAGNETIC SERVICE/capture certificate and therefore is not a
counterexample to eventual entry.

## 2. Source-uniform entry bound for the dominant AW block

The existing theorem contract does NOT currently prove a source-uniform
release bound small enough for AW entry into the r_FA ball.

The reason is quantitative and already independently audited. The previous
pointwise physical-AW tracking route is explicitly REFUTED on an admitted A21
history:

    sup ||a_hat_w-a_phys|| = 7.647 m/s^2 on a 16-s window.   (RE4)

By contrast, the radius-local FA proof at r=.15 uses

    ||a_hat_w-a_phys|| < 4.06(.15)=.609 m/s^2.               (RE5)

Therefore MARINE MOTION + IMU BIAS + MAGNETIC SERVICE do not imply the
pointwise .609-m/s^2 AW tube before local storage entry. The filter covariance
ceiling P_aw,aw<=16.48 I is an implication FROM storage to component error; it
does not bound deterministic AW mean error without a storage/action bound.
There is also no shipping hard projection of the AW mean analogous to the BA
or gyro-bias projections.

Consequently the desired source-uniform statement

    H18/release assumptions => V_o<=.15^2                    (RE6)

cannot be proved from the currently established capture/release lemmas. This
is not merely lack of a sharp constant: the required pointwise AW bridge used
to derive such an entry is known false on the admitted class.

The strongest existing source-uniform results relevant to AW are instead
signed/windowed and action based (Corollary A*, exact AW loop identities,
S-chain cancellation, complete chronological covariance/action identities).
They do not imply pointwise release localization.

## 3. Consequence for the proof architecture

The radius-local FA lemma remains correct once the local ball is reached, but
it cannot by itself bootstrap arbitrary H18/A21 release into that ball.
Attempting to use it for entry is circular.

The next controlling entry proof must therefore use one of two structures
already present in the repository, without adding assumptions:

1. a SHAPED outer region whose AW condition is the weaker signed/windowed
   quantity actually controlled by the shipping AW/S dynamics, followed by
   LaSalle/compactness to reach V<=.15^2; or
2. the global same-history FA reachability calculation
   (literal required physical acceleration -> private adaptation chronology ->
   linked AW correction/covariance action) to exclude the field-axis equality
   set before local storage entry.

A pointwise source-uniform AW tracking lemma should NOT be attempted again.

The immediate mathematical blocker is now precise: construct an entry
functional compatible with the proved signed/windowed AW control. Full
release V_o and pointwise AW error are not such functionals.
