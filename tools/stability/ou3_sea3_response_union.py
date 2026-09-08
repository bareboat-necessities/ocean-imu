"""SEA3+ response contract: unchanged linear vessels UNION correlated Stokes.

This is a model-domain extension, not a certificate that a captured word has
the required source root, Normal-Live bounds, bias history or finite storage.
Finite Stokes atoms are a new physical-model branch, never a quadrature proof
of the existing continuum linear-vessel branch. No global atom-count cap is
introduced: a Stokes source cell fixes its finite N and all root coordinates.
"""
from __future__ import annotations

from fractions import Fraction

LINEAR = "LINEAR_VESSEL"
STOKES = "STOKES_WAVE_FOLLOWING"
BRANCHES = (LINEAR, STOKES)
QUALIFICATION = "OU3_SEA3_PLUS_RESPONSE_UNION_V1"


def build() -> dict:
    return {
        "qualification": QUALIFICATION,
        "set_operation": "UNION",
        "branches": list(BRANCHES),
        "linear_vessel_domain": "tools/stability/ou3_sea3_directional_response_domain.json",
        "linear_vessel_domain_unchanged": True,
        "branch_selected_at_source_root_and_retained_at_every_prefix": True,
        "per_event_branch_switching_allowed": False,
        "legacy_continuum_replaced_by_finite_grid": False,
        "stokes": {
            "model": "superposed third-order bound harmonics; not an exact water-wave PDE solution",
            "atom_count": "any finite N, fixed on a source lineage; 128 is the pinned generator instance",
            "fundamental_frequency_hz": [0.02, 0.8],
            "order": 3,
            "component_steepness_upper_exact": "1/5",
            "coefficients": ["a", "k*a^2/2", "3*k^2*a^3/8"],
            "dispersion": "omega=2*pi*f; k=omega^2/g",
            "phase_graph": "theta_i(x,y,t)=k_i*d_i dot (x,y)-omega_i*t+phi_i; harmonic n uses n*theta_i",
            "directions_and_phases": "all root directions and phases, not seed 42 only",
            "amplitudes": "nonnegative common-root amplitudes; PM/JONSWAP quadrature and its Stokes normalization are source coordinates, not the continuum identity",
            "joint_output": "orbital displacement/velocity/acceleration and slopes at the displaced surface point; world-to-body attitude from those slopes",
            "gyro": "exact vee(-Cdot*C^T); the generator centered difference is a separately recorded measurement-discretization defect",
            "steady_drift": "retain constant Stokes velocity drift separately from centered orbital position; no fictitious drift acceleration",
            "hard_bound": "B_r=sum_i sum_n (n*omega_i)^r*c_in; ||d||<=B_0, ||v_orb||<=B_1, ||a||<=B_2",
            "slope_bound": "||s||<=sum_i sum_n n*k_i*c_in, also at the displaced point",
            "global_compact_primitive_cover_closed": False,
            "same_history_outward_output_cover_closed": False,
        },
        "shared_word_premises": [
            "unchanged Normal-Live acceleration/body-rate/chart bounds and vector PE",
            "same-history exact frontend/tuner/commit/scheduler and actual anisotropic R_S",
            "all shipping events and source-generated Live covariance",
            "BIAS0/1 true-bias root, forcing split and corrected-error recurrence",
        ],
        "linear_RAO_moment_and_surface_period_lemmas_scope": [LINEAR],
        "stokes_period_moments": "derive from the full phase-locked output; never substitute the fundamental surface Tz",
        "P3_extension_rule": "shared Normal-Live preconditions imply the same response-independent H18/A21 matrix certificate on each branch",
        "P4_extension_rule": "endpoint contraction AND finite every-prefix gain AND retention for both modes on EVERY branch/source leaf",
        "physical_deployment_inclusion_closed": False,
    }


def validate(d: dict) -> list[str]:
    # A changed union requires an explicit theorem/contract revision, not an
    # extra label silently accepted by an old certificate.
    return [] if d == build() else ["SEA3+ response union contract changed or incomplete"]


def stokes_hard_bounds(atoms) -> dict:
    """Exact rational algebra on supplied (omega,k,a) coordinates of ONE cell.

    Dispersion/phase provenance is a separate source premise. This helper does
    not claim transcendental validation of 2*pi*f or complete source admission.
    The point API audit below is explicitly floating point, not this theorem.
    """
    rows = [tuple(map(Fraction, row)) for row in atoms]
    if not rows:
        raise ValueError("Stokes cell must contain atoms (zero amplitudes are allowed)")
    bounds = [Fraction(0) for _ in range(3)]
    slope = Fraction(0)
    for omega, k, a in rows:
        if omega <= 0 or k <= 0 or a < 0 or k * a > Fraction(1, 5):
            raise ValueError("invalid Stokes atom or component steepness")
        coeffs = (a, k * a * a / 2, 3 * k * k * a**3 / 8)
        for n, c in enumerate(coeffs, 1):
            for r in range(3):
                bounds[r] += (n * omega)**r * c
            slope += n * k * c
    return {"B0_exact": str(bounds[0]), "B1_exact": str(bounds[1]),
            "B2_exact": str(bounds[2]), "slope_bound_exact": str(slope)}


def branch_mode_coverage(mode_checks: dict) -> dict:
    """Lift a response-independent implication; never accept partial modes."""
    if set(mode_checks) != {"H18", "A21"}:
        raise ValueError("both H18 and A21 must be covered")
    if any(type(v) is not bool for v in mode_checks.values()):
        raise ValueError("coverage needs Boolean certificate results")
    return {branch: dict(mode_checks) for branch in BRANCHES}


def all_branches_modes_closed(coverage: dict) -> bool:
    return (set(coverage) == set(BRANCHES)
            and all(set(coverage[b]) == {"H18", "A21"}
                    and all(v is True for v in coverage[b].values())
                    for b in BRANCHES))
