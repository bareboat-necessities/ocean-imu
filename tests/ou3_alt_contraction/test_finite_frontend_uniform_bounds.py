from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import unittest

from tools.stability.ou3_alt_contraction import finite_frontend_uniform_bounds as U
from tools.stability.ou3_alt_contraction import finite_band_binary32_contraction as BR
from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as FS
from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as S
from tools.stability.ou3_alt_contraction import finite_stats_binary32_runtime as ST


class Tests(unittest.TestCase):
    def test_all_time_induction_and_large_eager_variance_intermediate(self):
        doc = U.build()
        U.validate(doc)
        self.assertTrue(all(x > 0 for x in doc['strict_inductive_margins'].values()))
        self.assertFalse(doc['startup_deadline_required'])
        self.assertFalse(doc['startup_capture_and_complete_MEKF_totality_proved_here'])
        self.assertEqual(doc['sigma_target_upper'], 4)
        self.assertGreater(doc['largest_ordinary_intermediate_abs_upper'], F(10**19))
        self.assertLess(doc['largest_ordinary_intermediate_abs_upper'], F(2*10**19))
        doc['strict_inductive_margins']['p11'] = 1
        with self.assertRaisesRegex(ValueError, 'certificate changed'):
            U.validate(doc)

    def test_actual_compiled_configuration_and_mutated_setters(self):
        kwargs = dict(band_cfg=U.default_band_config(), stats_cfg=U.default_stats_config(),
            still_cfg=U.STILL_CFG.Config(), cutoff_hz=F(6), dt=U.DT,
            bench_noise_sigma=U.BENCH_SIGMA)
        U.require_config(**kwargs)
        for name, value in (('cutoff_hz', F(1)), ('dt', F(1, 200)),
                            ('bench_noise_sigma', F(0))):
            with self.assertRaises(ValueError):
                U.require_config(**(kwargs | {name: value}))
        with self.assertRaisesRegex(ValueError, 'shipping constants'):
            U.require_config(**(kwargs | {'stats_cfg': replace(kwargs['stats_cfg'], K_periods=8)}))

    def test_local_machine_covariance_all_contractions_at_extreme_box_faces(self):
        # Box-face checks audit the analytical rounding charge; the universal
        # induction itself uses exact polynomial inequalities, not this sample.
        for al, ah, sign in product(U.LOW_ALPHA, U.HIGH_ALPHA, (-1, 1)):
            predecessor = BR.State(U.LOWPASS, U.BAND, F(2), F(2*sign), F(8), True)
            out = BR.step(predecessor, x=U.VERTICAL, q_low=1-al, q_high=1-ah)
            for values, cap in ((out.lowpass_values, U.LOWPASS), (out.band_values, U.BAND),
                    (out.p00_values, F(2)), (out.p01_values, F(2)), (out.p11_values, F(8))):
                self.assertTrue(all(abs(x) < cap for x in values))

    def test_variance_readout_at_minimum_ready_weight_is_finite(self):
        weight = U.B.add(ST.READY_WEIGHT, F(1, 2**43))
        state = ST.State(mean_value=U.MEAN, mean_weight=weight,
                         sq_value=U.SQUARE, sq_weight=weight)
        self.assertTrue(state.var_ready)
        self.assertEqual(ST.variance_outcomes(state), (F(0),))
        self.assertTrue(U.B.is_binary32(U.B.mul(U.B.div(U.MEAN, weight),
                                               U.B.div(U.MEAN, weight))))

    def test_runtime_exp_profile_accepts_target_error_outside_ideal_RNE_cell(self):
        from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E
        from tools.stability.ou3_alt_contraction import finite_band_coefficients_binary32 as BC
        def shifted_exp(argument):
            lo, hi = E.exp_minus_enclosure(argument)
            ideal = U.B.rn32((lo+hi)/2)
            value = U.B.add(ideal, F(1, 2**24))
            self.assertFalse(E._interval_hits_rne_cell(lo, hi, value))
            self.assertTrue(E._interval_hits_exp_error_cell(lo, hi, value))
            return value
        f = U.B.rn32(F(1, 5))
        low = U.B.mul(F(1, 2), f); high = U.B.mul(F(4), f)
        arg = lambda x: U.B.mul(U.B.mul(BC.TWO_PI, x), U.DT)
        coeff = BC.produce(U.default_band_config(), f_ref=f, dt=U.DT,
            exp_low=shifted_exp(arg(low)), exp_high=shifted_exp(arg(high)))
        self.assertTrue(U.LOW_ALPHA[0] <= coeff.alpha_low <= U.LOW_ALPHA[1])
        stats = ST.coefficients(U.default_stats_config(), frequency=f, dt=U.DT,
                                exp_decay=shifted_exp(U.B.div(U.DT, F(20))))
        self.assertTrue(U.STATS_ALPHA[0] <= stats.alpha <= U.STATS_ALPHA[1])
        previous = VS.LPFState(F(1), F(6), True, 1)
        decay = shifted_exp(arg(F(6)))
        value = U.B.fma(decay, previous.value, U.B.mul(1-decay, F(2)))
        result = VS.lpf_step(previous, x=F(2), dt=U.DT, alpha_exp=decay, successor=value)
        self.assertEqual(result.state.value, value)

    def test_same_event_binder_and_carried_nonzero_state(self):
        import test_finite_machine_vertical_stillness_source as VHELP
        import test_finite_machine_frontend_sigma_source as FHELP
        initial = VS.begin(V.State(initialized=True), 0, False, S.State(), cutoff_hz=6)
        guarded = VHELP.quiet_guarded()
        vc, sc = VHELP.vertical_cfg(), VHELP.still_cfg()
        mah = VHELP.MAH.step_initialized(initial.vertical, vc, dt=U.DT,
            gyro=guarded.raw_gyro_body, acc=guarded.conditioned_accel_body)
        energy, attenuation = VHELP.still_operands(mah.vertical.vertical_accel, U.DT, sc)
        vertical = VS.step(initial, guarded, vc, sc, dt=U.DT,
            still_energy_successor=energy, still_attenuation_exp=attenuation)
        front = FHELP.source_step(FS.initial(), band_cfg=U.default_band_config(),
            stats_cfg=U.default_stats_config(), frequency=U.B.rn32(F(1, 5)),
            bench_noise_sigma=U.BENCH_SIGMA, x=vertical.band_input)
        U.require_event(vertical, front)
        self.assertGreater(front.state.stats.mean_weight, 0)
        self.assertGreater(front.state.band.machine.p11, 0)
        bad = replace(front.state, band=replace(front.state.band,
            machine=replace(front.state.band.machine, p11=9)))
        with self.assertRaisesRegex(ValueError, 'covariance outside'):
            U.require_state(vertical.state, bad)
        with self.assertRaises(TypeError):
            U.require_event(object(), front)


if __name__ == '__main__':
    unittest.main()
