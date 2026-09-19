"""Complement storage, conditional supply algebra and non-promotion tests."""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal as D, localcontext
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_theorem import h18_iss as iss
from tools.stability.ou3_theorem import h18_superword as base
from test_ou3_h18_superword import synthetic_export, rows, scaled_identity, STEP


def held_export(gains=(1.0, 0.5), bias_coupling=0.0):
    export = synthetic_export(samples=len(gains)-1, mag_prefix=len(gains)-1)
    maps = []
    for k, gain in enumerate(gains):
        a = scaled_identity(gain)
        for i in range(18, 21):
            a[i][i] = 1.0
        a[0][18] = k*bias_coupling
        maps.append(rows([[STEP*v for v in row] for row in a]))
    export['error_difference'] = maps
    export['root_difference_coarse'] = maps[0][:]
    export['endpoint_difference_coarse'] = maps[-1][:]
    return export


class SpectrumTests(unittest.TestCase):
    def test_rotated_symmetric_problem(self):
        with localcontext() as ctx:
            ctx.prec = 50
            value, residual = iss.symmetric_max([[D('.53'), D('.28')], [D('.28'), D('.53')]])
            self.assertLess(abs(value-D('.81')), D('1e-40'))
            self.assertLess(residual, D('1e-35'))

    def test_clustered_dominant_eigenvalues(self):
        with localcontext() as ctx:
            ctx.prec = 60
            value, _ = iss.symmetric_max([[D('0.99980000001'), D('1e-12')],
                                          [D('1e-12'), D('0.9998')]])
            self.assertGreater(value, D('0.99980000001'))
            self.assertLess(value, D('0.99980000002'))

    def test_bad_and_unconverged_spectra_refused(self):
        for matrix in ([[D('NaN')]], [[D(1), D(2)], [D(0), D(1)]], []):
            with self.subTest(matrix=matrix), self.assertRaises(base.SuperwordExportError):
                iss.symmetric_max(matrix)
        with self.assertRaises(base.SuperwordExportError):
            iss.symmetric_max([[D(1), D('.1')], [D('.1'), D(1)]], sweeps=0)


class SupplyTests(unittest.TestCase):
    def test_scalar_schur_gain_is_exact(self):
        with localcontext() as ctx:
            ctx.prec = 50
            gamma = iss.linear_bias_gain([[D('.5')]], [[D(2)]], D('.75'))
            self.assertLess(abs(gamma-D(6)), D('1e-40'))
            self.assertEqual(iss.young_multiplier(D('.25'), D('.75')), D('1.5'))

    def test_schur_gain_bounds_cross_term(self):
        with localcontext() as ctx:
            ctx.prec = 50
            m = [[D('.6'), D('.1')], [D(0), D('.5')]]
            b = [[D(1)], [D('.3')]]
            rho = D('.8')
            gamma = iss.linear_bias_gain(m, b, rho)
            for i in range(-3, 4):
                for j in range(-3, 4):
                    for u in (-2, 0, 2):
                        y = [D(i), D(j)]
                        z = [v+row[0]*u for v, row in zip(base.matvec(m, y), b)]
                        gap = base.dot(z, z) - rho*base.dot(y, y) - gamma*u*u
                        self.assertLessEqual(gap, D('1e-35'))

    def test_noncontractive_and_nonfinite_inputs_refused(self):
        for alpha, rho in (('1', '.9'), ('.9', '.8'), ('.5', '1'), ('NaN', '.8')):
            with self.subTest(alpha=alpha), self.assertRaises(base.SuperwordExportError):
                iss.young_multiplier(D(alpha), D(rho))
        with self.assertRaises(base.SuperwordExportError):
            iss.linear_bias_gain([[D(1)]], [[D(1)]], D('.9'))


class ComplementTests(unittest.TestCase):
    def test_identity_bias_is_an_input_not_a_contraction_coordinate(self):
        report = iss.evaluate(held_export(bias_coupling=.1), precision=40)
        self.assertAlmostEqual(report['endpoint_operator_ratio'], .25, places=12)
        self.assertGreater(report['linear_bias_supply_gain'], 0)
        self.assertEqual(report['prefix_operator_ratio_max'], 1)
        self.assertTrue(report['every_prefix_all_directions_evaluated'])

    def test_interior_prefix_expansion_not_hidden_by_endpoint(self):
        report = iss.evaluate(held_export((1, 1.5, .5)), precision=40)
        self.assertAlmostEqual(report['endpoint_operator_ratio'], .25)
        self.assertAlmostEqual(report['prefix_operator_ratio_max'], 2.25)
        self.assertEqual(report['prefix_of_maximum'], 1)

    def test_scale_discrepancy_can_exceed_an_observed_margin(self):
        export = held_export((1, .9999))
        coarse = export['endpoint_difference_coarse'][:]
        coarse[1] = STEP*.02
        export['endpoint_difference_coarse'] = coarse
        report = iss.evaluate(export, precision=40)
        self.assertLess(report['endpoint_operator_ratio'], 1)
        self.assertFalse(report['scale_discrepancy_fits_margin'])
        self.assertFalse(report['discrepancy_is_a_rigorous_error_bound'])

    def test_nonheld_or_coupled_exports_refused(self):
        export = held_export()
        export['error_difference'][-1][18*21] = .001
        with self.assertRaises(base.SuperwordExportError):
            iss.evaluate(export)
        export = held_export()
        export['covariance_upper'][0][18] = .001
        with self.assertRaises(base.SuperwordExportError):
            iss.evaluate(export)

    def test_no_outcome_promotes_finite_error_or_theorem(self):
        for gain in (.5, 1, 1.2):
            report = iss.evaluate(held_export((1, gain)), precision=40)
            for field in ('certificate_complete', 'source_uniform', 'obligation_discharged',
                          'finite_error_remainder_bounded', 'reference_forcing_bounded_in_storage'):
                self.assertFalse(report[field])
            if gain >= 1:
                self.assertIsNone(report['diagnostic_rho'])
                self.assertIsNone(report['linear_bias_supply_gain'])


class MagneticPhaseTests(unittest.TestCase):
    def test_missing_or_post_correction_phase_refused(self):
        for phase in (None, 'post_correction'):
            export = held_export()
            export['applied_magnetic_corrections'][0]['phase'] = phase
            with self.assertRaises(base.SuperwordExportError):
                base.validate(export)

    def test_nonfinite_pre_response_and_coarse_map_refused(self):
        export = held_export()
        export['applied_magnetic_corrections'][0]['pre_error_difference'][0] = float('nan')
        with self.assertRaises(base.SuperwordExportError):
            base.validate(export)
        export = held_export()
        export['endpoint_difference_coarse'][0] = float('inf')
        with self.assertRaises(base.SuperwordExportError):
            base.validate(export)

    def test_duplicate_fractional_and_unordered_event_prefixes_refused(self):
        for prefixes in ((1, 1), (2, 1), (1.5,)):
            export = held_export((1, .5, .25))
            event = export['applied_magnetic_corrections'][0]
            export['applied_magnetic_corrections'] = []
            for prefix in prefixes:
                copy = deepcopy(event)
                copy['prefix'] = prefix
                export['applied_magnetic_corrections'].append(copy)
            with self.assertRaises(base.SuperwordExportError):
                base.validate(export)

    def test_service_uses_pre_response_even_when_post_response_is_zero(self):
        export = held_export((1, .5, .25))
        event = export['applied_magnetic_corrections'][0]
        export['applied_magnetic_corrections'] = []
        for t in (1, 2):
            e = deepcopy(event)
            e['prefix'] = t
            a = scaled_identity(1)
            a[2][5] = t
            e['pre_error_difference'] = rows([[STEP*v for v in row] for row in a])
            export['applied_magnetic_corrections'].append(e)
        with localcontext() as ctx:
            ctx.prec = 40
            inv = base.inverse(base.square_from_rows(export['error_difference'][0], 21))
            zero = base.zeros(21, 21)
            maps = [base.SuperwordMap(zero, inv) for _ in range(3)]
            information, _ = base.applied_magnetic_information(export, maps, {
                'heading_coordinate_scale_rad': 1, 'axial_gyro_bias_coordinate_scale_rad_s': 1})
            self.assertGreater(information, 1)


if __name__ == '__main__':
    unittest.main()
