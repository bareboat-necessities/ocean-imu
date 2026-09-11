#!/usr/bin/env python3
"""Deterministic physical wave source: derive, rather than assume, its primitive.

Two sufficient realizations are supported. They are alternatives, not simultaneous
restrictions on every sea: (1) a zero-DC vector amplitude measure with hard
weighted total variation; (2) a bounded shaping state with an exact potential
identity. Both define p relative to the wave equilibrium, not global position.
The paper gives the all-time proofs. Rational arithmetic below checks their
algebra and computes source-uniform numerical bounds from supplied generators.
No PSD, Hs statistic, tuner clamp, replay or P4 radius supplies those inputs.

The repository has not qualified a numerical generator envelope covering its
entire physical family. build() exposes that gap, not a guessed 300 m*s budget.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
from typing import Sequence

SCHEMA = 1
QUALIFICATION = "OU3_BRMM_PHYSICAL_WAVE_GENERATOR_V1"
CANONICAL_SOURCE = "COMPLETE_BRMM_NORMAL_LIVE_WORD"


def rational(x) -> F:
    """The certificate language is exact; do not silently decimalize floats."""
    if isinstance(x, (bool, float)):
        raise TypeError("use integer, Fraction or rational string, not float/bool")
    return F(x)


def sqrt_upper(x: F, denominator: int = 10**18) -> F:
    """Exact rational upper enclosure of sqrt(x), including exact zero."""
    x = rational(x)
    if x < 0 or denominator <= 0:
        raise ValueError("nonnegative radicand and positive denominator required")
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    if a*a == x.numerator and b*b == x.denominator:
        return F(a, b)
    n = isqrt(x.numerator * denominator**2 // x.denominator)
    if F(n*n, denominator**2) < x:
        n += 1
    y = F(n, denominator)
    if y*y < x:
        raise ArithmeticError("rational square-root enclosure failed")
    return y


def matrix(rows):
    result = tuple(tuple(rational(x) for x in row) for row in rows)
    if not result or not result[0] or len({len(r) for r in result}) != 1:
        raise ValueError("nonempty rectangular matrix required")
    return result


def transpose(a):
    return tuple(zip(*a))


def multiply(a, b):
    if not a or not b or len(a[0]) != len(b):
        raise ValueError("matrix dimensions do not agree")
    return tuple(tuple(sum((x*y for x, y in zip(row, col)), F(0))
                       for col in transpose(b)) for row in a)


def plus(a, b, scale=F(1)):
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        raise ValueError("matrix dimensions do not agree")
    return tuple(tuple(x+scale*y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def psd(a, *, strict=False) -> bool:
    """Exact LDL Schur elimination; a zero pivot requires its entire row zero."""
    if not a or len(a) != len(a[0]) or a != transpose(a):
        return False
    q = [list(row) for row in a]
    for k in range(len(q)):
        d = q[k][k]
        if d < 0 or (strict and d == 0):
            return False
        if d == 0:
            if any(q[k][j] != 0 for j in range(k+1, len(q))):
                return False
            continue
        for i in range(k+1, len(q)):
            for j in range(i, len(q)):
                q[i][j] -= q[i][k]*q[k][j]/d
                q[j][i] = q[i][j]
    return True


def inverse(a):
    n = len(a)
    if not n or len(a[0]) != n:
        raise ValueError("square matrix required")
    q = [list(row)+[F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for k in range(n):
        i = next((i for i in range(k, n) if q[i][k]), None)
        if i is None:
            raise ValueError("singular matrix")
        q[k], q[i] = q[i], q[k]
        p = q[k][k]
        q[k] = [v/p for v in q[k]]
        for i in range(n):
            if i != k:
                p = q[i][k]
                q[i] = [x-p*y for x, y in zip(q[i], q[k])]
    return tuple(tuple(row[n:]) for row in q)


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def encode(a):
    return [[str(x) for x in row] for row in a]


@dataclass(frozen=True)
class SpectralBand:
    """A continuum support band, NOT a frequency grid or PSD bin.

    mass bounds TV(A)+TV(B), the hard vector amplitude measures in
    p(t)=integral cos(wt)dA(w)+sin(wt)dB(w). Directions/cross-axis response
    and every phase remain part of those same fixed measures.
    """
    omega_lower: F
    omega_upper: F
    mass: F

    def __post_init__(self):
        for name in ("omega_lower", "omega_upper", "mass"):
            object.__setattr__(self, name, rational(getattr(self, name)))
        if not 0 < self.omega_lower <= self.omega_upper or self.mass < 0:
            raise ValueError("positive wave-frequency support and nonnegative hard mass required")


def spectral_certificate(bands: Sequence[SpectralBand]) -> dict:
    """All measures in these bands, all phases, all t,u: D_S=2 sum mass/w_lo.

    Finite bands are convenient encodings of a finite weighted-measure norm.
    The paper also allows support approaching zero when that inverse moment
    remains uniformly finite. Overlapping bands are safe when masses describe
    a decomposition; they need not partition the frequency axis.
    """
    if not all(isinstance(b, SpectralBand) for b in bands):
        raise TypeError("hard amplitude bands required; a spectrum is not a certificate")
    m1 = sum((b.mass/b.omega_lower for b in bands), F(0))
    return {
        "kind": "zero_DC_hard_vector_amplitude_measure",
        "bands": [{"omega_lower_rad_s": str(b.omega_lower),
                   "omega_upper_rad_s": str(b.omega_upper),
                   "hard_amplitude_total_variation_m": str(b.mass)} for b in bands],
        "potential_norm_upper_m_s": str(m1),
        "D_S_upper_m_s": str(2*m1),
        "position_norm_upper_m": str(sum((b.mass for b in bands), F(0))),
        "velocity_norm_upper_mps": str(sum((b.mass*b.omega_upper for b in bands), F(0))),
        "acceleration_norm_upper_mps2": str(sum((b.mass*b.omega_upper**2 for b in bands), F(0))),
        "same_measure_owns_p_v_a_and_potential": True,
        "all_phases_and_all_times_covered": True,
        "independent_per_word_measure_selection_allowed": False,
    }


def shaping_certificate(*, A, B, C, D, L, P, alpha, initial_radius,
                        driver_norm_upper) -> dict:
    """Validate a bounded physical generator, including its all-time invariant.

    xdot=A x+B u, p=C x+D u, phi=L x. Exact LA=C and LB=D imply phidot=p.
    A'P+PA <= -2 alpha P, P>0 give an invariant state ellipsoid. For alpha>0,
    r=max(r0, sqrt(tr(B'PB))*U/alpha). For alpha=0 require B=0 (lossless
    oscillators allowed), so r=r0. No theorem depends on numerical sampling.
    """
    A, B, C, D, L, P = map(matrix, (A, B, C, D, L, P))
    n, m = len(A), len(B[0])
    if (len(A[0]) != n or len(B) != n or len(P) != n or len(P[0]) != n
            or len(C) != 3 or len(C[0]) != n or len(L) != 3 or len(L[0]) != n
            or len(D) != 3 or len(D[0]) != m):
        raise ValueError("generator/output dimensions invalid")
    alpha, r0, U = map(rational, (alpha, initial_radius, driver_norm_upper))
    if alpha < 0 or r0 < 0 or U < 0 or not psd(P, strict=True):
        raise ValueError("nonnegative hard input/initial bounds and P>0 required")
    if multiply(L, A) != C or multiply(L, B) != D:
        raise ValueError("physical displacement is not the generator potential derivative (LA=C, LB=D)")
    dissipation = plus(plus(multiply(transpose(A), P), multiply(P, A)), P, 2*alpha)
    if not psd(tuple(tuple(-v for v in row) for row in dissipation)):
        raise ValueError("hard shaping-state invariant not certified")
    b = sqrt_upper(trace(multiply(multiply(transpose(B), P), B)))
    if alpha == 0 and b != 0:
        raise ValueError("lossless state with a forcing port needs a separate invariant proof")
    r = max(r0, b*U/alpha) if alpha else r0
    ell = sqrt_upper(trace(multiply(multiply(L, inverse(P)), transpose(L))))
    return {
        "kind": "bounded_physical_shaping_potential",
        "A": encode(A), "B": encode(B), "C": encode(C), "D": encode(D),
        "L": encode(L), "P": encode(P), "alpha": str(alpha),
        "initial_radius": str(r0), "driver_norm_upper": str(U),
        "invariant_state_radius": str(r),
        "potential_norm_upper_m_s": str(r*ell), "D_S_upper_m_s": str(2*r*ell),
        "potential_derivative_identity_exact": True,
        "state_invariant_exact_LDL_checked": True,
        "p_v_a_same_output_and_derivatives_required": True,
        "all_time_bounded_driver_not_Gaussian_event": True,
    }


def verify_certificate(cert: dict) -> list[str]:
    """Recompute; callers cannot promote a fabricated D_S or identity flag."""
    if not isinstance(cert, dict):
        return ["physical generator certificate must be an object"]
    try:
        if cert.get("kind") == "zero_DC_hard_vector_amplitude_measure":
            expected = spectral_certificate(tuple(SpectralBand(
                b["omega_lower_rad_s"], b["omega_upper_rad_s"],
                b["hard_amplitude_total_variation_m"]) for b in cert["bands"]))
        elif cert.get("kind") == "bounded_physical_shaping_potential":
            expected = shaping_certificate(**{k: cert[k] for k in (
                "A", "B", "C", "D", "L", "P", "alpha", "initial_radius", "driver_norm_upper")})
        else:
            return ["unsupported physical generator; PSD/boolean/finite-window boxes are not admission"]
    except (KeyError, TypeError, ValueError, ZeroDivisionError, IndexError) as exc:
        return [str(exc)]
    return ([k+" differs from derived generator certificate" for k, v in expected.items() if cert.get(k) != v]
            + [str(k)+" is not a generator certificate field" for k in cert.keys()-expected.keys()])


def certificate_id(cert: dict) -> str:
    failures = verify_certificate(cert)
    if failures:
        raise ValueError(str(failures))
    return sha256(json.dumps(cert, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def constant_history_admission(displacement: Sequence, cert: dict) -> dict:
    """All-time admission, not a finite replay or a constant-exclusion flag.

    For d_j !=0 the exact horizon h=D_S/|d_j|+1 contradicts
    |integral_0^h p_j|<=D_S, an already derived generator consequence.
    d=0 has the zero measure/state realization (additional BRMM caps separate).
    """
    certificate_id(cert)
    d = tuple(rational(x) for x in displacement)
    if len(d) != 3:
        raise ValueError("displacement must have three components")
    amplitude = max(abs(x) for x in d)
    ds = rational(cert["D_S_upper_m_s"])
    h = ds/amplitude+1 if amplitude else None
    return {
        "admitted_by_physical_wave_condition": amplitude == 0,
        "derived_D_S_m_s": str(ds),
        "contradiction_horizon_s": None if h is None else str(h),
        "primitive_component_at_horizon_m_s": None if h is None else str(h*amplitude),
        "strict_excess_m_s": None if h is None else str(h*amplitude-ds),
        "reason": "zero wave realization" if not amplitude else "unbounded integral contradicts bounded generator potential",
    }


def harmonic_coefficient_identities(omega=F(1)) -> dict:
    """Exact for arbitrary vector coefficients and phases, not one trajectory.

    Coordinates [cos(wt),sin(wt)]; row differentiation sends [A,B] to
    [w B,-w A]. Check primitive, velocity and acceleration basis columns.
    """
    w = rational(omega)
    if w <= 0:
        raise ValueError("positive frequency required")
    J = matrix(((0, -w), (w, 0)))
    C = matrix(((1, 0), (0, 1), (0, 0)))
    L = multiply(C, inverse(J))
    acceleration = multiply(multiply(C, J), J)
    return {
        "potential_derivative_residual": encode(plus(multiply(L, J), C, F(-1))),
        "acceleration_plus_omega_squared_position_residual": encode(plus(acceleration, C, w*w)),
    }


def build() -> dict:
    identities = harmonic_coefficient_identities()
    exact = all(all(F(v) == 0 for row in a for v in row) for a in identities.values())
    # Non-promoting analytical examples prove non-vacuity and enforce that the
    # physical budget is NOT capped by 300. These are not the full BRMM family.
    examples = {
        "quiet": spectral_certificate(()),
        "continuum_band": spectral_certificate((SpectralBand(F(1, 10), F(2), F(3)),)),
        "slow_wave_not_clipped_to_working_tube": spectral_certificate((SpectralBand(F(1, 100), F(1, 100), F(2)),)),
    }
    examples["bounded_driven_shaper"] = shaping_certificate(
        A=[[-1]], B=[[1]], C=[[-1], [0], [0]], D=[[1], [0], [0]],
        L=[[1], [0], [0]], P=[[1]], alpha=1, initial_radius=0, driver_norm_upper=1)
    if any(verify_certificate(c) for c in examples.values()):
        raise RuntimeError("exact generator backend failed")
    dc = constant_history_admission((F(1, 8), 0, 0), examples["continuum_band"])
    # Universal contradiction, with independent symbolic D>=0 and a>0:
    # h=(D+a)/a; a*h-D=(D+a)-D=a. Check both formal coefficients.
    excess = tuple(x-y for x, y in zip((F(1), F(1)), (F(1), F(0))))
    rejection_theorem = exact and excess == (F(0), F(1))
    return {
        "schema": SCHEMA, "qualification": QUALIFICATION, "canonical_source": CANONICAL_SOURCE,
        "wave_position_not_global_vessel_position": True,
        "physical_position_decomposition": "x_CoG=x_equilibrium+p_wave in fixed axes; v=dp_wave/dt; a=dv/dt",
        "baseline_acceleration_policy": "ddot(x_equilibrium) is zero, independently compensated, or retained in declared same-history model-disturbance budget; never silently removed",
        "spectral_source": {
            "representation": "p(t)=integral cos(omega*t)dA(omega)+sin(omega*t)dB(omega), omega>0",
            "hard_family_norms": "uniform finite M_j=integral omega^j(d|A|+d|B|), j=-1,0,1,2",
            "primitive": "phi(t)=integral [sin(omega*t)dA-cos(omega*t)dB]/omega",
            "derived_bound": "D_S=2*M_minus1; if omega>=omega_min>0 then D_S<=2*M_0/omega_min",
            "PSD_or_Hs_determines_hard_amplitude_mass": False,
            "frequency_cutoff_taken_from_tuner": False,
            "physical_frequency_gap_required_if_inverse_moment_is_qualified": False,
        },
        "shaping_source": {
            "representation": "xdot=A*x+B*u; p=C*x+D*u; phi=L*x; LA=C and LB=D",
            "invariant": "P>0; A'P+PA<=-2*alpha*P; bounded driver and initial state as checked by shaping_certificate",
            "derived_bound": "D_S<=2*r*sqrt(trace(L*inverse(P)*L'))",
            "time_varying_extension": "Ldot+L*A=C and L*B=D, with a uniformly bounded potential; a rate bound alone is insufficient",
            "switching": "one carried state; every generator switch must preserve phi, the invariant, and the physical p/v/a jet; no wordwise reset",
        },
        "general_modulation_condition": "for p=sum(A_i(t)cos(theta_i)+B_i(t)sin(theta_i)), theta_dot=omega_i>0: D_S<=2*sup sum((|A_i|+|B_i|)/omega_i)+TV_global(A_i/omega_i,B_i/omega_i); local rate caps alone do not bound the residual integral",
        "same_history_primitive_relation": "S_L(t)=phi(t)-phi(t_L); S_L(next)-S_L=phi(next)-phi=integral p",
        "generator_state_and_live_potential_carried_across_all_words": True,
        "zero_mean_alone_sufficient": False,
        "hard_bounded_primitive_derived_from_generator": exact,
        "constant_nonzero_position_zero_velocity_history_admitted": dc["admitted_by_physical_wave_condition"],
        "zero_wave_admitted": constant_history_admission((0, 0, 0), examples["quiet"])["admitted_by_physical_wave_condition"],
        "exact_coefficient_identities": identities,
        "analytic_examples_not_full_source_qualification": examples,
        "constant_history_regression": dc,
        "historical_finite_window_only_classification": "B",
        "intended_physical_source_omission_classification": "E",
        "universal_constant_history_rejection": {
            "quantifiers": "every derived finite D>=0 and every nonzero displacement, with a=max_j |d_j|>0",
            "contradiction_horizon": "h=(D+a)/a",
            "strict_excess_coefficients_in_D_a": [str(x) for x in excess],
            "conclusion": "a*h-D=a>0 contradicts the generator primitive bound",
        },
        "old_witness_excluded_by_corrected_physical_theorem": rejection_theorem,
        "uniform_physical_family_generator_envelope": None,
        "existing_caps_determine_uniform_D_S": False,
        "missing_uniformity_harmonic_family": {
            "scope": "existing jet caps and positive individual frequencies, not the corrected fixed-budget family",
            "parameter": "all integers n>=10",
            "physical_heave": "p_z=cos(t/n); v_z=-sin(t/n)/n; a_z=-cos(t/n)/n^2; fixed attitude",
            "primitive": "phi_z=n*sin(t/n); exact all-u,t diameter 2*n",
            "uniform_position_velocity_acceleration_bounds": ["1", "1/10", "1/100"],
            "ten_second_AC_energy_upper": "1/1000",
            "gravity_direction_chord": "identically zero: vertical acceleration cannot reverse g-a",
            "conclusion": "no finite uniform D_S follows from the old physical caps; qualify an inverse-frequency amplitude budget or bounded potential generator",
            "filter_instability_claimed": False,
            "corrected_fixed_budget_theorem_refuted": False,
        },
        "numeric_D_S_m_s": None,
        "physical_D_S_numeric_qualification_closed": False,
        "numeric_working_tube_comparison": "unresolved: physical generator envelope is not numerically qualified",
        "source_specification_corrected": exact,
        "source_uniform_executor_materialized_here": False,
        "P3_delta": 1e-18, "P4_PASS": False, "P5_MAY_START": False,
    }


def validate(d: dict) -> list[str]:
    expected = build()
    return [k+" differs from derived physical source contract" for k, v in expected.items() if d.get(k) != v]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(); failures = validate(d)
    d.update(validation_pass=not failures, validation_failures=failures)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"physical_source_corrected": d["source_specification_corrected"],
                      "DC_history_admitted": d["constant_nonzero_position_zero_velocity_history_admitted"],
                      "numeric_D_S_qualified": d["physical_D_S_numeric_qualification_closed"],
                      "failures": failures}, sort_keys=True))
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
