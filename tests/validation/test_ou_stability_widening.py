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

    def test_analytical_foundation_has_explicit_vector_information_bound(self):
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

    def test_envelope_is_structural_not_a_second_quantitative_certificate(self):
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
        self.assertIn("not a parallel basin certificate", flat)
        for retired in (
            r"\label{eq:widen-kappaN}",
            r"\label{lem:widen-constructive-ues}",
            r"\label{eq:widen-Mrho}",
            r"\label{thm:widen-explicit-P}",
            r"\label{eq:widen-ell-star}",
        ):
            self.assertNotIn(retired, proof)
        self.assertNotIn("fallback", flat.casefold())

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
        flat = _flat(proof).casefold()
        self.assertIn("exact for arbitrary finite", flat)
        self.assertIn("no common-radius or componentwise radius relaxation", flat)
        self.assertIn("no global mekf claim", flat)

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

    def test_phase_e_uses_same_funnel_for_hybrid_and_model_mismatch(self):
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
        flat = _flat(proof)
        self.assertIn("same source-node Lyapunov metric", flat)
        self.assertNotIn(r"p_+^{\rm an}", proof)
        self.assertNotIn(r"\rho^{", proof)

    def test_phase_f_uses_product_manifold_funnel_drift(self):
        proof = _read(DOC / "w3d-stability-widening-phase-f.tex-part")
        for marker in (
            r"\label{eq:widen-stochastic-product-distance}",
            r"\label{eq:widen-stochastic-remainder}",
            r"\label{eq:widen-gaussian-moments}",
            r"\label{eq:widen-stochastic-margin}",
            r"\label{eq:widen-stochastic-floor}",
            r"\label{thm:widen-stochastic-drift}",
            r"\label{eq:widen-stochastic-coefficients}",
            r"\label{eq:widen-stochastic-exit}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("never subtracts two attitudes", flat)
        self.assertIn("same $W_i$, source graph, word family, and geodesic funnel", flat)
        self.assertNotIn(r"\kappa_N", proof)


if __name__ == "__main__":
    unittest.main()
