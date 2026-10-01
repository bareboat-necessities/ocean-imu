# Radius-local field-alignment exclusion (controlling LaSalle lemma)

This note is the controlling local exclusion for the MARINE dissipativity /
LaSalle route. It uses the literal zero-dissipation geometry and the retained
storage ball; it does not bound arbitrary innovations or invoke the scalar
vertical frontend.

## Lemma RL-A — storage-local nominal/physical acceleration relation

On every regular A21 prefix in the candidate set V=e'P^-1 e <= r^2 for which
the inherited AW marginal ceiling P_aw,aw <= 16 I is applicable, matrix
Cauchy--Schwarz gives

    ||e_aw|| = ||a_hat_w-a_phys|| <= 4 r.                    (RL-A1)

The inherited active-BA marginal P_ba,ba <= I/1600 likewise gives

    ||e_ba|| <= r/40.                                        (RL-A2)

These are component consequences of the SAME full covariance/storage,
including all cross covariance. They are not pointwise AW-tracking physical
assumptions.

For the field-alignment lemma below RL-A2 is not added to RL-A1. The literal
accelerometer attitude Jacobian is formed from the nominal CoG force
R_wb(a_hat_w-g); BA is a separate measurement column and the lever-arm term is
attitude-independent in this Jacobian. On the zero-dissipation set the
measurement residual itself is zero. Adding BA, lever and fast sensor boxes to
RL-A1 would therefore double-count terms that are absent from the geometric
condition.

## Lemma RL-B — sampled physical force/field separation

Let b be the applicable fixed unit committed magnetic reference and
P_B=I-bb'. On the qualified field-fraction branch

    sigma_w=||e_z x b|| >= 1/5,

so with the shipping gravity magnitude g=9.80665 m/s^2,

    ||P_B g|| >= g/5 = 1.96133 m/s^2.                       (RL-B1)

For one MARINE physical history, v'=a and ||v||<=Vmax=5.5 m/s. Define

    d(t)=P_B(a(t)-g).

On any interval W=[t,t+T],

    integral_W d(s) ds
      = P_B[v(t+T)-v(t)] - T P_B g.

Hence

    max_(s in W) ||d(s)||
      >= ||P_B g|| - 2 Vmax/T
      >= 1.96133 - 11/T.                                    (RL-B2)

Physical acceleration is 100-Lipschitz. Gravity and P_B are fixed in this
committed-reference calculation, hence d is also 100-Lipschitz. On the
retained regular A21 branch actually applied accelerometer epochs have gap at
most h_acc=.006 s. Choose a maximizer s*. An applied epoch k exists within
h_acc, so

    ||d(t_k)||
      >= 1.96133 - 11/T - 100(.006)
      = 1.36133 - 11/T
      =: m_phys(T).                                         (RL-B3)

This deliberately uses J*h_acc, not the sharper trapezoidal J*h_acc/4 charge:
RL-B3 transfers one continuous maximizer to one actually applied epoch.

Why "actually applied" is legitimate on the retained mathematical branch:
the wrapper invokes the accelerometer correction at every valid MEKF-driven
IMU sample; R_acc has a strictly positive floor, so the innovation covariance
is SPD in the real-arithmetic retained domain. Invalid/nonfinite/failed-solve
branches are arithmetic/domain obligations and are not silently counted as
applied service.

RL-B3 is positive for

    T > 11/1.36133 = 8.080... s.                            (RL-B4)

## Theorem RL-FA — radius-local exclusion of persistent field alignment

Suppose a nonzero forward-complete zero-dissipation MARINE trajectory remains
inside V<=r^2. The existing zero-dissipation classification reduces its only
nonzero candidate to the field-axis attitude mode. Persistence requires at
every applied accelerometer epoch

    P_B(a_hat_w-g)=0.                                       (RL-FA1)

Using RL-A1,

    ||P_B(a_phys-g)||
      = ||P_B(a_phys-a_hat_w)||
      <= 4r.                                                (RL-FA2)

But RL-B says every T-window contains an actually applied accelerometer epoch
with

    ||P_B(a_phys-g)|| >= m_phys(T).                          (RL-FA3)

Therefore persistence is impossible whenever

    M(r,T):=1.96133 - 11/T - .6 - 4r > 0.                  (RL-FA4)

At the already used 17-s nuisance/root warm-up,

    m_phys(17)
      = 1.96133 - 11/17 - .6
      = 0.714271176470588... m/s^2,

so every

    r < r_FA,max(17)
      = m_phys(17)/4
      = 0.178567794117647...                                (RL-FA5)

is admissible for this exclusion. A deliberately conservative explicit choice
is

    r_FA = 0.15,     T_FA = 17 s,                           (RL-FA6)

for which

    M(.15,17)=0.114271176470588... m/s^2 > 0.              (RL-FA7)

The six-degree local attitude domain is a separate coordinate-domain
restriction; r_FA is the dimensionless covariance-storage radius. RL-A1
already supplies the physical units for 4r.

### Defect bookkeeping

No epsilon_B, lever, calibration, fast-accelerometer, attitude, or BA term is
added to RL-FA4:

* the explicit committed-reference premise sigma_w>=1/5 is already the
  magnetic geometry used in RL-B1; using g cos80 instead is an alternative,
  weaker branch, not an additional simultaneous defect;
* the zero-dissipation candidate is stated in the nominal WORLD CoG geometry
  P_B(a_hat_w-g)=0, so no finite attitude-coordinate conversion is needed;
* lever arm is attitude-independent in the implemented CoG attitude Jacobian;
* BA is a separate measurement column and its r/40 bound is not part of
  a_hat_w-a_phys;
* calibration/fast sensor residuals affect the measurement record, but zero
  dissipation makes the corresponding innovation/action zero; they do not
  change the nominal Jacobian vector whose collinearity defines RL-FA1.

If a later theorem weakens RL-FA1 to approximate/nonzero-dissipation
collinearity, those terms must be restored through the literal dissipation
inequality. They must not be inserted into the exact invariant-set lemma.

Consequently, on the stated retained regular MARINE branch,

    Inv_MARINE({D=0}) intersect {V<=.15^2} = {0},            (RL-FA8)

conditional only on the already-proved zero-dissipation classification and
the retained-domain component/cadence premises cited above.

## Compactness consequence

Let K be the normalized compact retained MARINE history/root class
V_root=1, V_prefix<=r_FA^2 after homogeneous scaling, with the closed event
strata and same-history continuation already used by the variational proof.
D>=0 is continuous/lower-semicontinuous on each closed literal stratum.

If no finite m and eta_D>0 existed with

    sum_(j=0)^(m-1) D_(k+j) >= eta_D V_k,                   (RL-C1)

then for every n there would be a normalized admitted n-window with cumulative
dissipation tending to zero. Compactness/diagonal extraction gives a
forward-complete admitted limit with D=0 at every operation. RL-FA8 forces
that normalized limit to be zero, contradicting V_root=1. Hence some finite
m and eta_D>0 exist. Storage telescoping then gives

    V_(k+m) <= (1-eta_D) V_k                                (RL-C2)

for the homogeneous retained system.

RL-C1 is an existence result; this argument does not manufacture a numerical
eta_D or close nonlinear/source supply, capture, transitions or float32
totality. Those remain subsequent obligations.
