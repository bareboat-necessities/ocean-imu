from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"tools/stability"))
import ou3_brmm_physical_wave_source as W
import ou3_brmm_wave_primitive_binding as BIND
import ou3_brmm_centered_primitive_transition as PRIM


def shaper(**changes):
    # xdot=-x+u, p=xdot, phi=x. This is an actual bounded-input high-pass
    # displacement realization, not an S cap. Smooth output-jet caps remain
    # separate physical admission requirements.
    args = dict(A=[[-1]], B=[[1]], C=[[-1], [0], [0]], D=[[1], [0], [0]],
                L=[[1], [0], [0]], P=[[1]], alpha=1,
                initial_radius=0, driver_norm_upper=1)
    args.update(changes)
    return W.shaping_certificate(**args)


class PhysicalWaveSourceTests(unittest.TestCase):
    def test_positive_individual_frequencies_are_not_a_uniform_family_bound(self):
        # For every n>=10 the same p/v/a caps hold; the exact primitive diameter
        # is 2*n. This is a parameterized algebraic family, not sampled evidence
        # of failure of the corrected fixed-budget source or of the filter.
        for n in (10, 200, 10**50):
            cert = W.spectral_certificate((W.SpectralBand(F(1, n), F(1, n), 1),))
            self.assertEqual(F(cert["D_S_upper_m_s"]), 2*n)
            self.assertLessEqual(F(cert["velocity_norm_upper_mps"]), F(1, 10))
            self.assertLessEqual(F(cert["acceleration_norm_upper_mps2"]), F(1, 100))
        self.assertFalse(W.build()["existing_caps_determine_uniform_D_S"])

    def test_forged_live_potential_cannot_reintroduce_an_independent_S_ball(self):
        from dataclasses import replace
        cert = shaper()
        b = BIND.enter_live(certificate=cert, source_id="s", origin_id="L",
                            state_id="x", generator_state=(0,), p=(0,0,0), v=(0,0,0))
        # Even forging both S and the origin consistently must fail: the
        # original physical generator cannot have reached this Live potential.
        b = replace(b, live_potential=(F(1000), F(0), F(0)),
                    primitive=replace(b.primitive, S_L=(F(-1000), F(0), F(0))))
        with self.assertRaisesRegex(ValueError, "derived generator image"):
            BIND.validate_boundary(b, cert)

    def test_exact_sqrt_encloses_without_sampling(self):
        for x in (F(0), F(4, 9), F(2), F(1, 10**70), F(10**50+1, 3)):
            y = W.sqrt_upper(x)
            self.assertGreaterEqual(y*y, x)
        self.assertEqual(W.sqrt_upper(F(4, 9)), F(2, 3))
        with self.assertRaises(TypeError):
            W.sqrt_upper(0.5)

    def test_continuum_band_bound_is_source_uniform_exact(self):
        c = W.spectral_certificate((W.SpectralBand(F(1, 10), 2, 3),
                                     W.SpectralBand(2, 3, F(1, 2))))
        self.assertEqual(c["D_S_upper_m_s"], "121/2")
        self.assertEqual(c["velocity_norm_upper_mps"], "15/2")
        self.assertEqual(c["acceleration_norm_upper_mps2"], "33/2")
        self.assertTrue(c["all_phases_and_all_times_covered"])
        self.assertEqual(W.verify_certificate(c), [])
        for w in (0, -1):
            with self.assertRaises(ValueError):
                W.SpectralBand(w, 2, 3)

    def test_every_formal_harmonic_coefficient_cancels_exactly(self):
        for w in (F(1, 1000000), F(7, 11), F(100000)):
            d = W.harmonic_coefficient_identities(w)
            self.assertTrue(all(all(F(x) == 0 for row in v for x in row) for v in d.values()))

    def test_DC_witness_rejected_by_derived_contradiction_not_cap(self):
        for cert in (W.spectral_certificate((W.SpectralBand(1, 2, 3),)), shaper()):
            for sign in (-1, 1):
                d = W.constant_history_admission((F(sign, 8), 0, 0), cert)
                self.assertFalse(d["admitted_by_physical_wave_condition"])
                self.assertGreater(F(d["primitive_component_at_horizon_m_s"]), F(d["derived_D_S_m_s"]))
            self.assertTrue(W.constant_history_admission((0, 0, 0), cert)["admitted_by_physical_wave_condition"])
        d = W.build()
        self.assertFalse(d["constant_nonzero_position_zero_velocity_history_admitted"])
        self.assertTrue(d["old_witness_excluded_by_corrected_physical_theorem"])
        self.assertFalse(d["physical_D_S_numeric_qualification_closed"])

    def test_300_is_not_a_source_admission_parameter(self):
        c = W.spectral_certificate((W.SpectralBand(F(1, 100), F(1, 100), 2),))
        self.assertEqual(c["D_S_upper_m_s"], "400")
        self.assertEqual(W.verify_certificate(c), [])
        c["D_S_upper_m_s"] = "300"
        self.assertTrue(W.verify_certificate(c))

    def test_PSD_and_declaration_flags_do_not_qualify_a_generator(self):
        for c in ({"kind": "PSD", "Hs": 8.5},
                  {"kind": "boolean", "S_bounded": True},
                  {"kind": "finite_window", "position_bounded": True}):
            self.assertTrue(W.verify_certificate(c))
        cert = shaper(); cert["unmodeled_DC_position"] = "1"
        self.assertTrue(W.verify_certificate(cert))
        for cert in (None, [], {"kind": "bounded_physical_shaping_potential", "A": [[]]}):
            self.assertTrue(W.verify_certificate(cert))
        d = W.build(); d["physical_D_S_numeric_qualification_closed"] = True
        self.assertTrue(W.validate(d))

    def test_shaper_state_invariant_and_potential_are_derived(self):
        c = shaper()
        self.assertEqual(c["invariant_state_radius"], "1")
        self.assertEqual(c["D_S_upper_m_s"], "2")
        self.assertEqual(W.verify_certificate(c), [])
        with self.assertRaisesRegex(ValueError, "potential derivative"):
            shaper(D=[[0], [0], [0]])  # ordinary low-pass can have a DC output
        with self.assertRaisesRegex(ValueError, "invariant"):
            shaper(alpha=2)  # forged decay rate
        with self.assertRaisesRegex(ValueError, "lossless"):
            shaper(alpha=0)  # cannot claim boundedness from a lossless driven state
        c["L"][0][0] = "2"
        self.assertTrue(W.verify_certificate(c))

    def test_lossless_harmonic_state_covers_all_initial_phases(self):
        c = W.shaping_certificate(A=[[0, -1], [1, 0]], B=[[0], [0]],
                                 C=[[2, 0], [0, 0], [0, 0]], D=[[0], [0], [0]],
                                 L=[[0, 2], [0, 0], [0, 0]], P=[[1, 0], [0, 1]],
                                 alpha=0, initial_radius=1, driver_norm_upper=0)
        self.assertEqual(c["D_S_upper_m_s"], "4")
        self.assertEqual(W.verify_certificate(c), [])

    def test_bounded_local_modulation_rate_is_not_sufficient(self):
        # (1+eps*cos(t))*cos(t) = cos(t)+eps/2+(eps/2)*cos(2t).
        # Both amplitude and its rate are bounded, carrier frequency is 1,
        # but the integral has a nonzero linear term eps*t/2.
        eps = F(1, 10)
        dc_coefficient = eps/2
        self.assertEqual(dc_coefficient, F(1, 20))
        cert = W.spectral_certificate((W.SpectralBand(1, 2, 2),))
        self.assertFalse(W.constant_history_admission((dc_coefficient, 0, 0), cert)["admitted_by_physical_wave_condition"])
        self.assertFalse(W.build()["zero_mean_alone_sufficient"])

    def test_same_potential_is_carried_across_successive_words(self):
        c = shaper()
        b = BIND.enter_live(certificate=c, source_id="physical-1", origin_id="Live-1",
                            state_id="x0", generator_state=(0,), p=(0, 0, 0), v=(1, 0, 0))
        zero = (F(0), F(0), F(0))
        # One smooth local physical jet p(t)=t, a=0, phi(t)=t^2/2.
        # x=phi is generated by u=t+t^2/2 in the bounded shaper on [0,1/2].
        for k in range(1, 3):
            t = F(k, 4)
            b = BIND.advance(b, c, PRIM.MomentWitness("moment-"+str(k), zero, zero, zero),
                             next_generator_state=(t*t/2,), next_state_id="x"+str(k), h=F(1, 4))
            self.assertEqual(b.primitive.S_L[0], t*t/2)
            self.assertEqual(b.primitive.p[0], t)
            self.assertEqual(b.live_potential, zero)
            self.assertEqual(b.primitive.centered_S_origin_witness_id, "Live-1")
        with self.assertRaisesRegex(ValueError, "same history"):
            BIND.advance(b, c, PRIM.MomentWitness("m3", zero, zero, zero),
                         next_generator_state=(0,), next_state_id="x3", h=F(1, 4))
        with self.assertRaisesRegex(ValueError, "invariant"):
            BIND.advance(b, c, PRIM.MomentWitness("m3", zero, zero, zero),
                         next_generator_state=(2,), next_state_id="x3", h=F(1, 4))

    def test_finite_DC_prefix_not_excluded_by_fiat_but_cannot_continue_forever(self):
        c = shaper()
        b = BIND.enter_live(certificate=c, source_id="physical", origin_id="Live",
                            state_id="x0", generator_state=(0,), p=(F(1, 8), 0, 0), v=(0, 0, 0))
        zero = (F(0), F(0), F(0))
        for k in range(1, 8):
            b = BIND.advance(b, c, PRIM.MomentWitness(str(k), zero, zero, zero),
                             next_generator_state=(F(k, 8),), next_state_id=str(k), h=F(1))
        with self.assertRaisesRegex(ValueError, "invariant"):
            BIND.advance(b, c, PRIM.MomentWitness("later", zero, zero, zero),
                         next_generator_state=(F(9, 8),), next_state_id="later", h=F(2))
        # Endpoints alone are an outer relation, not complete source admission.
        # The exact all-time generator theorem supplies the independent rejection.
        self.assertFalse(W.constant_history_admission((F(1, 8), 0, 0), c)["admitted_by_physical_wave_condition"])

    def test_exact_LDL_rejects_indefinite_zero_pivot(self):
        self.assertFalse(W.psd(W.matrix([[0, 1], [1, 0]])))
        self.assertTrue(W.psd(W.matrix([[0, 0], [0, 1]])))
        self.assertFalse(W.psd(W.matrix([[0, 0], [0, 1]]), strict=True))


if __name__ == "__main__":
    unittest.main()
