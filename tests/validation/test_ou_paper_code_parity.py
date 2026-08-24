#!/usr/bin/env python3
"""Pin the deployed OU-II/OU-III paper claims to the implementation.

This is intentionally a source-level publication contract.  Numerical evidence
has its own provenance/evidence tests; these checks prevent the papers from
silently drifting away from the filter state layouts, selected adaptation laws,
constants, cadence semantics, and startup source path.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]


def text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


class OUPaperCodeParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ou3_core = text("src/kalman_ou_iii/Kalman3D_Wave_OU_III.h")
        cls.ou3_wrap = text("src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h")
        cls.ou2_core = text("src/kalman_ou_ii/Kalman3D_Wave_OU_II.h")
        cls.ou2_wrap = text("src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h")
        cls.ou3_state = text("doc/kalman_ou_iii/w3d-state.tex-part")
        cls.ou3_proc = text("doc/kalman_ou_iii/w3d-proc-cont.tex-part")
        cls.ou3_meas = text("doc/kalman_ou_iii/w3d-meas.tex-part")
        cls.ou3_adapt = text("doc/kalman_ou_iii/w3d-adaptation-motivation.tex-part")
        cls.ou3_impl = text("doc/kalman_ou_iii/w3d-fus-methods.tex-part")
        cls.ou3_init = text("doc/kalman_ou_iii/w3d-init.tex-part")
        cls.ou2_paper = text("doc/kalman_ou_ii/ou2-dual-regularization-mse.tex")

    def test_ou3_state_and_process_match_21_state_core(self):
        self.assertIn("static constexpr int EXT_ADD = 12", self.ou3_core)
        self.assertIn("static constexpr int OFF_S   = BASE_N + 6", self.ou3_core)
        self.assertIn("static constexpr int OFF_AW  = BASE_N + 9", self.ou3_core)
        self.assertIn("BASE_N + 12", self.ou3_core)

        self.assertIn("N_X=21", self.ou3_state)
        for token in (
            r"\delta\vct{\theta}",
            r"\delta\vct{b}_g",
            r"\delta\vct{v}",
            r"\delta\vct{p}",
            r"\delta\vct{S}",
            r"\delta\vct{a}_w",
            r"\delta\vct{b}_a",
        ):
            self.assertIn(token, self.ou3_state)
        self.assertIn(r"\dot{\vct{v}}=\vct{a}_w", self.ou3_proc)
        self.assertIn(r"\dot{\vct{p}}=\vct{v}", self.ou3_proc)
        self.assertIn(r"\dot{\vct{S}}=\vct{p}", self.ou3_proc)
        self.assertIn(r"\vct{z}_S=0=\vct{S}+\vct{n}_S", self.ou3_meas)

    def test_ou3_spectral_mse_default_matches_paper(self):
        self.assertIn(
            "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;",
            self.ou3_wrap,
        )
        self.assertRegex(
            self.ou3_wrap,
            r"R_S_MSE_COEFF_DEFAULT\s*=\s*0\.0538f",
        )
        self.assertIn(
            "0.0148f * 0.0148f * FREQ_SMOOTHER_DT",
            self.ou3_wrap,
        )
        self.assertIn("std::pow(u, 6.0f / 7.0f)", self.ou3_wrap)
        self.assertIn("/ std::sqrt(TS)", self.ou3_wrap)
        self.assertIn(
            "if (rs_law_ != RSAdaptationLaw::Cubic) return 1.0f;",
            self.ou3_wrap,
        )

        self.assertIn(r"C_J=0.0538", self.ou3_impl)
        self.assertIn(r"\widehat\sigma_{a,B,\star}^{\,6/7}", self.ou3_impl)
        self.assertIn(r"\tau_\star^{24/7}T_{S,\star}^{-1/2}", self.ou3_impl)
        self.assertIn(r"\sqrt{R_a}=\SI{0.0148}{m.s^{-2}}", self.ou3_adapt)
        self.assertIn("No additional reference-cadence factor is applied", self.ou3_impl)

    def test_ou3_cadence_and_bias_dynamics_match_paper(self):
        for token in (
            "PSEUDO_UPDATE_PERIOD_NOMINAL_S = 0.015f",
            "PSEUDO_UPDATE_TAU_NOMINAL_S = 1.1f",
            "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = FREQ_SMOOTHER_DT",
            "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.25f",
            "ADAPT_TAU_SEA_PERIODS          = 0.40f",
            "ADAPT_RS_MULT              = 1.5f",
        ):
            self.assertIn(token, self.ou3_wrap)
        self.assertIn(r"c_T=\frac{\SI{15}{ms}}{\SI{1.1}{s}}", self.ou3_impl)
        self.assertIn(r"T_{S,\min}=\SI{5}{ms}", self.ou3_impl)
        self.assertIn(r"T_{S,\max}=\SI{250}{ms}", self.ou3_impl)
        self.assertIn(r"c_{\mathrm{adapt}}=0.40", self.ou3_impl)
        self.assertIn(r"m_{r_S}=1.5", self.ou3_impl)
        self.assertIn(r"\tau_b=\SI{5000}{s}", self.ou3_impl)

    def test_ou3_magnetic_refinement_uses_independent_proxy_tilt(self):
        self.assertIn(
            "const Eigen::Quaternionf q_tilt_bw = impl_.startupProxyTiltQuat();",
            self.ou3_wrap,
        )
        self.assertIn("Magnetic acquisition uses proxy tilt rather than the MEKF", self.ou3_init)
        self.assertIn("later\nindependent acquisition refines that reference", self.ou3_init)

    def test_ou2_state_and_dual_pseudos_match_core(self):
        self.assertIn("static constexpr int EXT_ADD = 9", self.ou2_core)
        self.assertIn("static constexpr int OFF_P   = BASE_N + 3", self.ou2_core)
        self.assertIn("static constexpr int OFF_AW  = BASE_N + 6", self.ou2_core)
        self.assertIn("measurement_update_position_pseudo", self.ou2_core)
        self.assertIn("measurement_update_velocity_pseudo", self.ou2_core)

        self.assertIn(
            "local error coordinates containing attitude, gyroscope bias, velocity,\n"
            "displacement, latent OU acceleration, and accelerometer bias",
            self.ou2_paper,
        )
        self.assertIn("one on world-frame\ndisplacement $p$ and one on world-frame velocity $v$", self.ou2_paper)
        self.assertIn("both applied with period $T_S$", self.ou2_paper)

    def test_ou2_physical_mse_default_matches_paper(self):
        self.assertIn(
            "PseudoAdaptationLaw pseudo_law_ = PseudoAdaptationLaw::PhysicalMSE;",
            self.ou2_wrap,
        )
        self.assertRegex(
            self.ou2_wrap,
            r"R_PSEUDO_MSE_COEFF_DEFAULT\s*=\s*0\.1116f",
        )
        self.assertRegex(
            self.ou2_wrap,
            r"R_PSEUDO_MSE_RATIO_DEFAULT\s*=\s*0\.4611f",
        )
        self.assertIn(
            "ACC_NOISE_FLOOR_SIGMA_DEFAULT * ACC_NOISE_FLOOR_SIGMA_DEFAULT * FREQ_SMOOTHER_DT",
            self.ou2_wrap,
        )
        self.assertIn("std::pow(u, 0.8f) / std::sqrt(TS)", self.ou2_wrap)
        self.assertIn("r_v = r_p / (ratio * tau);", self.ou2_wrap)
        self.assertIn(
            "if (pseudo_law_ != PseudoAdaptationLaw::Empirical) return 1.0f;",
            self.ou2_wrap,
        )

        for token in (
            r"C_P=0.1116",
            r"\frac{C_P}{C_V}=0.4611",
            r"\sigma_{\rm floor}=\SI{0.12}{m.s^{-2}}",
            r"h=\SI{0.005}{s}",
            r"\tau^{12/5}T_S^{-1/2}",
            r"\tau^{7/5}T_S^{-1/2}",
        ):
            self.assertIn(token, self.ou2_paper)
        self.assertIn(
            "no additional historical\ncadence renormalization is applied",
            self.ou2_paper,
        )

    def test_ou2_residual_noise_parameter_is_not_a_live_tuner_alias(self):
        # The two defaults share one characterization, but the implementation
        # deliberately stores them as separate knobs.  The paper must say so.
        self.assertIn(
            "float acc_noise_floor_sigma_ = ACC_NOISE_FLOOR_SIGMA_DEFAULT;",
            self.ou2_wrap,
        )
        self.assertIn(
            "float pseudo_accel_noise_density_ = R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT;",
            self.ou2_wrap,
        )
        self.assertIn("setPseudoAccelNoiseDensity", self.ou2_wrap)
        self.assertIn("stored as a separate scheduler parameter", self.ou2_paper)
        self.assertIn("does not live-read the tuner's noise-floor setting", self.ou2_paper)


if __name__ == "__main__":
    unittest.main()
