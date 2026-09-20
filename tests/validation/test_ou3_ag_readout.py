from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.ag_readout import (
    bootstrap, certificate, coefficient_relaxation_obstruction, exact_readout,
    factor_rows, noise_action_lower, readout_action, supplied_fixture,
)
from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import (
    add, identity, is_psd, matmul, transpose,
)


def posterior(p, events):
    """Independent exact Joseph reference; only measurement inverses (<=3)."""
    for event in events:
        if event['kind'] in ('prediction', 'reset'):
            f = event['F'] if event['kind'] == 'prediction' else event['G']
            p = matmul(matmul(f, p), transpose(f))
            if event['kind'] == 'prediction':
                u = event['U']
                p = add(p, matmul(u, transpose(u)))
        elif event.get('accepted', True):
            h, v = event['H'], event['V']
            r = matmul(v, transpose(v))
            s = add(matmul(matmul(h, p), transpose(h)), r)
            k = matmul(matmul(p, transpose(h)), inverse(s))
            a = add(identity(len(p)), matmul(k, h), F(-1))
            p = add(matmul(matmul(a, p), transpose(a)), matmul(matmul(k, r), transpose(k)))
    return p


class HistoricalReadoutTests(unittest.TestCase):
    def test_exact_21_state_dominance_with_unbounded_AG_prior_and_cross_terms(self):
        events = supplied_fixture()
        action = readout_action(events, exact_readout(events), identity(15))['action']
        for root_scale in (1, 10**6):
            c = identity(21)
            for i in range(6):
                c[i][i] = F(root_scale)
                c[i][6+i] = F(1, 4)
            p = matmul(c, transpose(c))
            end = posterior(p, events)
            self.assertTrue(is_psd(add(action, [row[:6] for row in end[:6]], F(-1))))
            self.assertTrue(any(any(row[6:]) for row in end[:6]))

    def test_uncancelled_root_cannot_be_hidden_as_small_numeric_residual(self):
        events = supplied_fixture()
        reader = exact_readout(events)
        reader[0][0] += F(1, 10**60)
        with self.assertRaisesRegex(ValueError, 'uncancelled AG root'):
            readout_action(events, reader, identity(15))

    def test_fixed_reader_is_not_uniform_under_varying_coefficients(self):
        events = supplied_fixture()
        reader = exact_readout(events)
        events[4]['F'][0][3] += F(1, 10**6)
        with self.assertRaisesRegex(ValueError, 'uncancelled AG root'):
            readout_action(events, reader, identity(15))
        # A reader depending on the realized coefficients restores exact
        # cancellation; uniform bounded action still needs a separate proof.
        result = readout_action(events, exact_readout(events), identity(15))
        self.assertTrue(result['AG_root_cancelled_exactly'])
        self.assertFalse(result['source_uniform_verified'])

    def test_rejected_measurements_do_not_repair_rank(self):
        events = supplied_fixture()
        for index in (1, 5):
            events[index]['accepted'] = False
        # The two magnetic samples alone leave at least one AG direction.
        with self.assertRaisesRegex(ValueError, 'lacks full column rank'):
            exact_readout(events)

    def test_matrix_bootstrap_and_failed_margin(self):
        events = supplied_fixture()
        action = readout_action(events, exact_readout(events), identity(15))['action']
        p = events[0]
        result = bootstrap(action, identity(15), p['F'], p['U'], F(1, 10**10))
        self.assertTrue(is_psd(result['AG_loss_lower']))
        self.assertLess(result['rho_upper'], 1)
        self.assertFalse(result['source_uniform_verified'])
        with self.assertRaisesRegex(ValueError, 'not PSD'):
            bootstrap(action, identity(15), p['F'], p['U'], F(1, 100))
        with self.assertRaises(ValueError):
            bootstrap(action, identity(15), p['F'], p['U'], 0)

    def test_refuse_omitted_hard_events_and_singular_resets(self):
        events = supplied_fixture()
        for event in ({'kind': 'release'}, {'kind': 'reset', 'G': [[0]*21 for _ in range(21)]}):
            with self.assertRaises(ValueError):
                exact_readout(events+[event])

    def test_nuisance_cross_correlation_is_not_replaced_by_diagonal(self):
        events = supplied_fixture()
        upper = identity(15)
        # The factor reader uses the y-acceleration row in this fixture.
        upper[10][13] = upper[13][10] = F(1, 2)
        full = readout_action(events, exact_readout(events), upper)['action']
        diagonal = readout_action(events, exact_readout(events), identity(15))['action']
        self.assertNotEqual(full, diagonal)
        self.assertTrue(is_psd(full))

    def test_factor_pivot_does_not_invert_a_vanishing_first_component(self):
        for small in (F(0), F(1, 10**60)):
            rows = [[small, 0], [1, 0], [0, 1]]
            self.assertEqual(factor_rows(rows), [1, 2])

    def test_all_row_noise_lower_is_necessary_but_not_an_action_upper_bound(self):
        events = supplied_fixture()
        action = readout_action(events, exact_readout(events), identity(15))['action']
        lower = noise_action_lower(events)
        self.assertTrue(is_psd(add(action, lower, F(-1))))
        self.assertFalse(is_psd(add(lower, action, F(-1))))

    def test_exact_noise_enclosure_preserves_correlated_columns(self):
        from tools.stability.ou3_theorem.ag_readout_source_diagnostic import rational_upper_factor
        q = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        u = rational_upper_factor(q)
        self.assertTrue(u[1][0])
        enclosed = matmul(u, transpose(u))
        self.assertTrue(is_psd(add(enclosed, q, F(-1))))
        with self.assertRaises(ValueError):
            rational_upper_factor([[1, 2], [2, 1]])
        singular = [[F(1), F(1)], [F(1), F(1)]]
        factor = rational_upper_factor(singular)
        self.assertEqual(matmul(factor, transpose(factor)), singular)

    def test_coefficient_relaxation_is_singular_but_not_a_reachable_counterexample(self):
        report = coefficient_relaxation_obstruction()
        self.assertLess(F(report['nominal_aw_norm_squared']), F('8.8')**2)
        self.assertEqual(report['full_AG_array_rank'], 4)
        self.assertEqual(report['null_Rayleigh_margin'], '-1')
        self.assertFalse(report['LO_equals_terminal_AG_map_possible'])
        self.assertFalse(report['nominal_mean_recursion_satisfied'])
        self.assertFalse(report['shipping_counterexample'])

    def test_unreachable_prior_family_does_not_refute_shipping(self):
        report = certificate()
        obstruction = report['prior_scale_obstruction']
        self.assertEqual(F(obstruction['Rayleigh_margin_ceiling']), F(1, 10**12)-F(1, 10**6))
        self.assertFalse(obstruction['shipping_reachability_claimed'])
        self.assertFalse(obstruction['shipping_stability_refuted'])
        for key in ('uniform_historical_readout_verified', 'six_column_source_uniform_loss_verified',
                    'full_21_covariance_upper_verified', 'source_uniform_A21_linear_dissipativity',
                    'theorem_closed'):
            self.assertFalse(report[key])


if __name__ == '__main__':
    unittest.main()
