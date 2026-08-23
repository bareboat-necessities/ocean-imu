"""Publication/source contract for the OU--III analytical stability widening."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src" / "kalman_ou_iii"
COMMON = REPO_ROOT / "src" / "kalman_ou_common"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


class OUIIIStabilityWideningContractTests(unittest.TestCase):
    def test_all_phases_are_wired_into_manuscript(self):
        main = _read(DOC / "kalman_ou-w3d.tex")
        for include in (
            r"\input{w3d-analytical-stability-widening.tex-part}",
            r"\input{w3d-stability-widening-phase-c.tex-part}",
            r"\input{w3d-stability-widening-phase-d.tex-part}",
            r"\input{w3d-stability-widening-phase-e.tex-part}",
            r"\input{w3d-stability-widening-phase-f.tex-part}",
        ):
            self.assertIn(include, main)

    def test_phase_a_has_explicit_vector_information_bound(self):
        proof = _read(DOC / "w3d-analytical-stability-widening.tex-part")
        for marker in (
            r"\label{eq:widen-vector-packet}",
            r"\label{eq:widen-vector-mu}",
            r"\label{lem:widen-vector-coercivity}",
            r"\label{eq:widen-gamma-lower}",
            r"\label{lem:widen-two-packet-bg}",
            r"\label{eq:widen-alpha6}",
        ):
            self.assertIn(marker, proof)
        self.assertIn(r"1-\sqrt{1-s_{fm}^2}", proof)
        self.assertIn(r"\alpha_6^{\rm an}", proof)

    def test_phase_a_envelope_does_not_depend_on_adaptation_law(self):
        proof = _read(DOC / "w3d-analytical-stability-widening.tex-part")
        flat = _flat(proof)
        for marker in (
            r"\Pi^{\rm env}",
            r"\label{eq:widen-envelope-bounds}",
            r"\label{lem:widen-envelope-translation}",
            r"\label{thm:widen-envelope-ues}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("No EMA recurrence", flat)
        self.assertIn("does not depend on the adaptation-law exponents", flat)
        self.assertIn("performance mechanism rather than a stability-critical feedback law", flat)

    def test_phase_b_constructs_ues_and_lyapunov_constants(self):
        proof = _read(DOC / "w3d-analytical-stability-widening.tex-part")
        for marker in (
            r"\label{eq:widen-kappaN}",
            r"\label{lem:widen-constructive-ues}",
            r"\label{eq:widen-Mrho}",
            r"\label{thm:widen-explicit-P}",
            r"\label{eq:widen-P-bounds}",
            r"\label{eq:widen-P-decrease}",
            r"\label{eq:widen-pq-constants}",
            r"\label{eq:widen-ell-star}",
        ):
            self.assertIn(marker, proof)

    def test_phase_c_uses_exact_large_angle_group_dissipation(self):
        proof = _read(DOC / "w3d-stability-widening-phase-c.tex-part")
        for marker in (
            r"\label{eq:widen-SO3-energy}",
            r"\label{eq:widen-exact-source-correction}",
            r"\label{eq:widen-exact-injection}",
            r"\label{eq:widen-exact-energy-change}",
            r"\label{eq:widen-reset-G-bound}",
            r"\label{eq:widen-large-angle-sector}",
            r"\label{thm:widen-large-angle-dissipation}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("exact for arbitrary finite", flat)
        self.assertIn("no common-radius or componentwise radius relaxation", flat)
        self.assertIn("no global MEKF claim", flat)

        src = _read(SRC / "Kalman3D_Wave_OU_III.h")
        common = _read(COMMON / "KalmanOUCoreMath.h")
        self.assertIn("quat_from_delta_theta(dtheta)", src)
        self.assertIn("Identity() + T(0.5)*skew(dtheta)", common)

    def test_phase_d_reset_discards_prior_attitude_but_not_topology(self):
        proof = _read(DOC / "w3d-stability-widening-phase-d.tex-part")
        for marker in (
            r"\label{lem:widen-global-reset}",
            r"\label{eq:widen-global-reset-bound}",
            r"\label{eq:widen-reset-image}",
            r"\label{thm:widen-source-global-startup}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("for every pre-reset attitude", flat)
        self.assertIn("does not assert global asymptotic stability", flat)
        self.assertIn("no absolute-yaw claim", flat)

        src = _read(SRC / "Kalman3D_Wave_OU_III.h")
        self.assertIn("qref = quaternion_from_acc(acc_n);", src)

    def test_phase_e_closes_hybrid_and_model_mismatch_routes_analytically(self):
        proof = _read(DOC / "w3d-stability-widening-phase-e.tex-part")
        for marker in (
            r"\label{eq:widen-reset-to-set}",
            r"\label{eq:widen-cooldown-factor}",
            r"\label{eq:widen-hybrid-mu}",
            r"\label{thm:widen-hybrid-small-gain}",
            r"\label{eq:widen-model-mismatch}",
            r"\label{eq:widen-deterministic-drift}",
            r"\label{thm:widen-model-mismatch-iss}",
            r"\label{eq:widen-tau-mismatch-bound}",
        ):
            self.assertIn(marker, proof)

    def test_phase_f_derives_gate_safe_stochastic_drift(self):
        proof = _read(DOC / "w3d-stability-widening-phase-f.tex-part")
        for marker in (
            r"\label{eq:widen-stochastic-remainder}",
            r"\label{eq:widen-gaussian-moments}",
            r"\label{eq:widen-stochastic-margin}",
            r"\label{eq:widen-stochastic-floor}",
            r"\label{thm:widen-stochastic-drift}",
            r"\label{eq:widen-stochastic-coefficients}",
            r"\label{eq:widen-stochastic-exit}",
        ):
            self.assertIn(marker, proof)


if __name__ == "__main__":
    unittest.main()
