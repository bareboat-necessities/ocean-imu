"""Structural contract for the OU--III Phase-3 Live-handoff certificate.

The Phase-2 version of this file pinned prose and numbers and was removed for
being premature.  This one deliberately pins only what a rewrite of the section
must not silently lose: that the certificate is wired into the article and into
CI, that the argument is carried in the covariance metric rather than converted
out of it, that the four kinds of number the paper distinguishes are separate
objects in the code, and that a timeout-forced handoff cannot be reported as
theorem-certified.  Numbers and wording are left free.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc" / "kalman_ou_iii"
SRC = ROOT / "src"
TEST = ROOT / "tests" / "kalman_ou_iii"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


class OUIIILiveHandoffCertificateContractTests(unittest.TestCase):
    def test_article_wires_the_certificate_parts_in_order(self):
        # Phase 1 rides inside the semiglobal part; the three certificate
        # subsections are wired at the top level, and the constructive handoff
        # section has to come after the analytical closure it now rests on.
        semiglobal = _read(DOC / "w3d-semiglobal-stability.tex-part")
        self.assertIn(r"\input{w3d-block-local-iss.tex-part}", semiglobal)

        text = _read(DOC / "kalman_ou-w3d.tex")
        phase2 = r"\input{w3d-live-basin-certificate.tex-part}"
        analytic = r"\input{w3d-computer-assisted-live-basin.tex-part}"
        phase3 = r"\input{w3d-live-handoff-certificate.tex-part}"
        for marker in (phase2, analytic, phase3):
            self.assertIn(marker, text)
        self.assertLess(text.index(phase2), text.index(analytic))
        self.assertLess(text.index(analytic), text.index(phase3))

    def test_each_certificate_part_is_input_exactly_once(self):
        # Both this phase and the analytical closure wire the Phase-2 part;
        # inputting it twice would duplicate its labels and its text.
        joined = "\n".join(
            _read(p)
            for p in sorted(DOC.glob("*.tex-part")) + sorted(DOC.glob("*.tex"))
        )
        for part in (
            "w3d-live-basin-certificate.tex-part",
            "w3d-computer-assisted-live-basin.tex-part",
            "w3d-live-handoff-certificate.tex-part",
        ):
            self.assertEqual(
                joined.count(r"\input{" + part + "}"), 1, f"{part} is not input once"
            )

    def test_phase3_rests_on_the_analytical_closure(self):
        text = _read(DOC / "w3d-live-handoff-certificate.tex-part")
        self.assertIn(r"\ref{sec:iss-computer-assisted-certificate}", text)
        self.assertIn(r"\label{sec:cert-sources}", text)
        self.assertIn(r"\label{eq:cert-analytic-map}", text)

    def test_phase3_carries_the_argument_in_the_covariance_metric(self):
        text = _read(DOC / "w3d-live-handoff-certificate.tex-part")
        for marker in (
            r"\label{lem:cert-monotone}",       # one-sample metric gain <= 1
            r"\label{lem:cert-scale-invariance}",
            r"\label{thm:cert-metric-ues}",
            r"\label{thm:cert-metric-entrance}",
            r"\label{eq:cert-S-zero}",          # the integration-epoch identity
            r"\label{eq:cert-ceff}",
        ):
            self.assertIn(marker, text)

    def test_phase3_separates_proof_from_declaration(self):
        text = _read(DOC / "w3d-live-handoff-certificate.tex-part")
        self.assertIn("configured physical envelope", text)
        self.assertIn("exact by construction", text)
        self.assertIn(r"\label{sec:cert-handoff-bounds}", text)
        # The section has to state the gap rather than only the improvement.
        self.assertIn(r"\label{sec:cert-numbers}", text)

    def test_certificate_distinguishes_proved_from_measured_constants(self):
        src = _read(SRC / "kalman_ou_iii" / "LiveEntranceCertificate.h")
        self.assertIn("enum class CertificateSource", src)
        self.assertIn("Analytical", src)
        self.assertIn("IntervalMeasured", src)
        self.assertIn("intervalMeasured()", src)
        self.assertIn("theoremCertified", src)
        # The default has to be the proved set, or a passing certificate would
        # silently mean the weaker thing.
        self.assertIn("CertificateSource source = CertificateSource::Analytical;", src)

        wrapper = _read(SRC / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h")
        # Only a theorem-certified handoff may reach the Certified state.
        self.assertIn("live_certificate_.theoremCertified()", wrapper)
        self.assertIn("passed-unproved-constants", wrapper)

    def test_certificate_header_keeps_the_four_kinds_of_number_apart(self):
        src = _read(SRC / "kalman_ou_iii" / "LiveEntranceCertificate.h")
        for marker in (
            "enum class BoundSource",
            "Exact,",
            "Envelope,",
            "Measured,",
            "External",
            "struct LiveEnvelope",          # declared, not measured
            "struct LiveBasinConstants",    # interval constants
            "struct LiveHandoffObservables",
            "weakestSource",
            "LiveCertFailure",
        ):
            self.assertIn(marker, src)

    def test_a_failed_certificate_is_readable(self):
        src = _read(SRC / "kalman_ou_iii" / "LiveEntranceCertificate.h")
        for field in (
            "basin_lhs",
            "basin_rhs",
            "margin",
            "nonlinear_radius",
            "small_gain_nu",
            "sensitive_metric_norm",
            "kinematic_metric_norm",
            "eta_gravity_rad",
            "eta_heading_rad",
        ):
            self.assertIn(field, src)

    def test_timeout_forced_live_cannot_be_reported_as_certified(self):
        src = _read(SRC / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h")
        self.assertIn("enum class LiveCertification", src)
        for state in ("NotLive", "Uncertified", "Certified"):
            self.assertIn(state, src)
        self.assertIn("liveCertification()", src)
        # Certification must not survive a Live-interval boundary.
        self.assertIn("liveIntervalEpoch", src)
        self.assertIn("live_cert_epoch_", src)

    def test_integration_epoch_reset_exists_on_the_filter(self):
        src = _read(SRC / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h")
        self.assertIn("reset_integral_epoch", src)
        self.assertIn("seed_world_accel", src)
        self.assertIn("seed_gyro_bias", src)
        self.assertIn("block_sigma", src)

    def test_diagnostic_gates_on_the_metric_constants(self):
        src = _read(TEST / "live_basin_diagnostic.cpp")
        # Fixed-size Eigen matrices must not ask for thin U/V.
        self.assertNotIn("ComputeThinU", src)
        self.assertNotIn("ComputeThinV", src)
        self.assertIn("all_reference_metric_contract", src)
        self.assertIn("alpha_max", src)
        self.assertIn("PHASE3_SUMMARY", src)

    def test_permanent_tests_are_in_the_normal_ci_path(self):
        run = _read(TEST / "run_tests.sh")
        makefile = _read(TEST / "Makefile")
        for target in ("live_basin_diagnostic",
                       "live_entrance_certificate-test",
                       "live_handoff_validation"):
            self.assertIn("./" + target, run)
            self.assertIn(target, makefile)


if __name__ == "__main__":
    unittest.main()
