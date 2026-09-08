"""Check physical generator admission independently of a replay's SEA3 label.

The exact rational inequality concerns the generating response model. It does
not exclude a different continuum source with the same finite output samples.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/stability"))
import ou3_sea3_response_union as UNION  # noqa: E402
MANIFEST = Path(__file__).with_name("fixtures") / "ou3_physical_generator_provenance.json"


def particle_response_test(frequency, gain, corner, power):
    f, g, fc, p = map(Fraction, (frequency, gain, corner, power))
    if not (f > 0 and g >= 0 and fc > 0 and p.denominator == 1 and p >= 2):
        raise ValueError("invalid response envelope")
    allowed_squared = (g * min(Fraction(1), (fc / f) ** int(p))) ** 2
    return {
        "frequency_hz_exact": str(f),
        "particle_response_norm_squared_exact": "2",
        "maximum_allowed_norm_squared_exact": str(allowed_squared),
        "squared_gap_exact": str(Fraction(2) - allowed_squared),
        "response_envelope_satisfied": Fraction(2) <= allowed_squared,
    }


def stokes_response_test(observed: dict) -> dict:
    """Pinned API point check, not outward certification of a complete word."""
    atoms = observed.get("atoms", [])
    if len(atoms) != observed.get("frequency_count") or not atoms:
        raise ValueError("missing full Stokes amplitude/frequency root")
    if observed.get("order") != 3:
        raise ValueError("Stokes order changed")
    gravity = observed["gravity_mps2"]
    if not math.isfinite(gravity) or gravity <= 0:
        raise ValueError("invalid gravity")
    bounds = [0.0, 0.0, 0.0]
    steepness = []
    for atom in atoms:
        f, a = atom["frequency_hz"], atom["amplitude_m"]
        if not (math.isfinite(f) and math.isfinite(a) and 0.02 - 1e-14 <= f <= 0.8 + 1e-14 and a >= 0):
            raise ValueError("Stokes atom outside declared fundamental domain")
        omega = 2 * math.pi * f
        k = omega**2 / gravity
        steepness.append(k * a)
        if k * a > 0.2:
            raise ValueError("Stokes steepness premise failed")
        c = (a, k*a*a/2, 3*k*k*a**3/8)
        for key, expected in zip(("second_coefficient_m", "third_coefficient_m"), c[1:]):
            if not math.isclose(atom[key], expected, rel_tol=2e-13, abs_tol=1e-25):
                raise ValueError("independent harmonic coefficient breaks Stokes root graph")
        for n, cn in enumerate(c, 1):
            for r in range(3):
                bounds[r] += (n * omega)**r * cn
    return {
        "response_model_admission": "ADMITTED",
        "branch": UNION.STOKES,
        "max_component_steepness_point": max(steepness),
        "hard_bound_formula_point_evaluation": dict(zip(("B0", "B1", "B2"), bounds)),
        "outward_complete_word_certificate": False,
        "source_code_provenance_required_for_phase_and_joint_output": True,
        "full_Normal_Live_word_membership": "UNDETERMINED",
    }


def build(generator_root: Path, observed: dict):
    manifest = json.loads(MANIFEST.read_text())
    commit = subprocess.check_output(
        ["git", "-C", str(generator_root), "rev-parse", "HEAD"], text=True).strip()
    if commit != manifest["generator_commit"]:
        raise ValueError("physical generator commit differs from provenance")
    for name, digest in manifest["generator_file_sha256"].items():
        if hashlib.sha256((generator_root / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"physical generator source changed: {name}")
    if "SIM_DATA_VERSION ?= v1.1.3" not in (ROOT / "Makefile").read_text():
        raise ValueError("deployed simulation data version changed")
    if observed["frequency_count"] != 128 or observed["order"] != 3:
        raise ValueError("observed generator family changed")
    if not (abs(observed["highest_fundamental_hz"] - 0.8) < 1e-12
            and observed["highest_third_harmonic_amplitude_m"] > 0):
        raise ValueError("expected nonzero highest harmonic is absent")
    response = json.loads((ROOT / "tools/stability/ou3_sea3_directional_response_domain.json").read_text())["response_contract"]
    # This corner dominates the unchanged LINEAR branch, not the new union.
    gain = str(response["peak_translation_gain_range"][1])
    corner = str(response["rolloff_corner_hz_range"][1])
    power = str(response["high_frequency_rolloff_power_min"])
    inequality = particle_response_test("12/5", gain, corner, power)
    stokes = stokes_response_test(observed)
    return {
        "qualification": "PINNED_PHYSICAL_GENERATOR_SEA3_ADMISSION_AUDIT",
        "provenance": manifest,
        "SEA3_response_union": UNION.build(),
        "Stokes_response_test": stokes,
        "public_generator_API_observation": observed,
        "exact_response_test": inequality,
        "fundamental_input_obstruction": "a linear response cannot create the nonzero 3*f_max harmonic above the generator's fundamental support",
        "full_Stokes_elevation_input_obstruction": "the surface-particle translational gain violates the declared SEA3 envelope at 12/5 Hz",
        "generator_has_one_correlated_phase_direction_history": True,
        "generator_is_complete_continuum_SEA3_provider": False,
        "linear_vessel_response_admission": "REJECTED" if not inequality["response_envelope_satisfied"] else "REQUIRES_REMAINING_PREMISES",
        "direct_generator_response_admission": stokes["response_model_admission"],
        "admission_scope": "RESPONSE_MODEL_ONLY_NOT_COMPLETE_RETAINED_WORD",
        "gyro_centered_difference_defect_certified": False,
        "finite_sample_membership_under_another_realization": "UNDETERMINED",
        "retained_payload_rebound_to_regenerated_history": False,
        "true_bias_noise_model_in_simulator": "turn-on offset plus random walk and white measurement noise; not an unforced OU history",
        "bias_mismatch_silently_charged_to_ISS": False,
        "canonical_P4_falsified": False,
        "P4_promoted": False,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--generator-root", type=Path, required=True)
    ap.add_argument("--observation", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    report = build(args.generator_root, json.loads(args.observation.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("PHYSICAL_SEA3_GENERATOR_ADMISSION", report["direct_generator_response_admission"],
          json.dumps(report["exact_response_test"], sort_keys=True))
    print("FINITE_SAMPLE_MEMBERSHIP", report["finite_sample_membership_under_another_realization"])


if __name__ == "__main__":
    main()
