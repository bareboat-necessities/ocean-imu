"""Semantic/source contract for the mode-aware OU--III Live proof."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src" / "kalman_ou_iii"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


class OUIIIStabilityPhaseAContractTests(unittest.TestCase):
    def test_live_proof_has_source_modes(self):
        proof = _read(DOC / "w3d-iss-stability.tex-part")
        flat = _flat(proof)
        for marker in (
            r"\label{eq:iss-held-state}",
            r"\label{eq:iss-active-state}",
            r"\label{eq:iss-held-ues}",
            r"\label{eq:iss-held-forced}",
            r"\label{eq:iss-21-ues}",
            r"\label{thm:iss-21-ues}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("active 18-state", flat)
        self.assertIn("active 21-state", flat)
        self.assertIn("frozen accelerometer-bias error", flat)
        self.assertNotIn("homogeneous full-heading 21-state linearized", flat)

    def test_source_reachable_schedule_is_not_cartesian(self):
        proof = _read(DOC / "w3d-iss-stability.tex-part")
        flat = _flat(proof)
        for marker in (
            r"\label{eq:iss-source-ema}",
            r"\label{eq:iss-source-slew}",
            r"\Pi^{\rm src}_{k,m}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("source-reachable family", flat)
        self.assertIn("not an arbitrary Cartesian product", flat)

    def test_deployed_aw_process_floor_closes_translational_process_route(self):
        proof = _read(DOC / "w3d-iss-stability.tex-part")
        src = _read(SRC / "SeaStateFusionFilter_OU_III.h")
        for marker in (
            r"\label{eq:iss-aw-process-floor}",
            r"\label{eq:iss-trans-FG}",
            r"\label{eq:iss-trans-controllability}",
            r"\label{eq:iss-trans-process-gramian}",
            r"\label{eq:iss-source-Q-lower}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("std::max(0.05f, band_noise_floor_sigma_())", src)
        self.assertIn(r"\det[", proof)
        self.assertIn("=-1", proof)

    def test_s_observability_credits_spread_firings(self):
        proof = _read(DOC / "w3d-iss-stability.tex-part")
        for marker in (
            r"\label{eq:iss-S-spread-index}",
            r"q_W\Delta_{S,-}",
            r"\label{eq:iss-S-smin}",
            r"\label{lem:iss-axis-uco}",
        ):
            self.assertIn(marker, proof)
        self.assertIn(r"(q_W\Delta_{S,-})^6", proof)

    def test_three_s_route_is_explicitly_conditional(self):
        proof = _read(DOC / "w3d-iss-stability.tex-part")
        flat = _flat(proof)
        for marker in (
            r"\label{eq:iss-three-S-gramian}",
            r"\label{lem:iss-three-s-coupled}",
            r"\mu_a>0",
        ):
            self.assertIn(marker, proof)
        self.assertIn("conditional vector-information bound is verified", flat)

    def test_bias_hold_matches_source(self):
        src = _read(SRC / "Kalman3D_Wave_OU_III.h")
        self.assertIn(
            "const T phi_b = acc_bias_updates_enabled_ ? std::exp(-Ts / tau_b) : T(1);",
            src,
        )
        self.assertIn("if (!use_ba) freeze_acc_bias_rows_(PCt);", src)
        self.assertIn("if (!use_ba) freeze_acc_bias_rows_(K);", src)

    def test_no_mag_stays_on_yaw_quotient(self):
        proof = _read(DOC / "w3d-iss-stability.tex-part")
        for marker in (
            r"\label{sec:iss-nomag-quotient}",
            r"\label{eq:iss-nomag-pe}",
            r"\label{thm:iss-nomag-quotient}",
            "No convergence of absolute yaw is claimed",
        ):
            self.assertIn(marker, proof)


if __name__ == "__main__":
    unittest.main()
