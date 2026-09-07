"""Check physical generator admission independently of a replay's SEA3 label.

The exact rational inequality concerns the generating response model. It does
not exclude a different continuum source with the same finite output samples.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
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
    # This corner dominates every response in the current declared family.
    gain = str(response["peak_translation_gain_range"][1])
    corner = str(response["rolloff_corner_hz_range"][1])
    power = str(response["high_frequency_rolloff_power_min"])
    inequality = particle_response_test("12/5", gain, corner, power)
    return {
        "qualification": "PINNED_PHYSICAL_GENERATOR_SEA3_ADMISSION_AUDIT",
        "provenance": manifest,
        "public_generator_API_observation": observed,
        "exact_response_test": inequality,
        "fundamental_input_obstruction": "a linear response cannot create the nonzero 3*f_max harmonic above the generator's fundamental support",
        "full_Stokes_elevation_input_obstruction": "the surface-particle translational gain violates the declared SEA3 envelope at 12/5 Hz",
        "generator_has_one_correlated_phase_direction_history": True,
        "generator_is_complete_continuum_SEA3_provider": False,
        "direct_generator_response_admission": "REJECTED" if not inequality["response_envelope_satisfied"] else "REQUIRES_REMAINING_PREMISES",
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
