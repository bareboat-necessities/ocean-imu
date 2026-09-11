"""Probe construction and fail-closed scope, without compiling in unit discovery."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_alt_contraction import shipping_finite_identity as S
from tools.stability.ou3_alt_contraction import proof_plan as PLAN


class ShippingFiniteIdentityTests(unittest.TestCase):
    def test_observation_erasure_recovers_exact_shipping_tokens(self):
        original = (ROOT/S.CORE).read_text()
        transformed, patches = S.instrument_core(original)
        recovered = S.erase_observation_calls(transformed)
        self.assertEqual(''.join(recovered.split()), ''.join(original.split()))
        self.assertGreater(len(patches), 10)
        self.assertEqual((ROOT/S.CORE).read_text(), original)
        self.assertLess(transformed.index('ou3_alt_probe::record(2, *this)'),
                        transformed.index('    apply_pending_aw_covariance_inflation_();',
                                          transformed.index('ou3_alt_probe::record(2, *this)')))

    def test_unknown_source_shape_fails_instead_of_guessing_hook(self):
        original = (ROOT/S.CORE).read_text()
        with self.assertRaises(ValueError):
            S.instrument_core(original.replace('    qref.normalize();', '    qref.normalize( );'))
        with self.assertRaises(ValueError):
            S.instrument_core(S.instrument_core(original)[0])

    def test_braces_in_comments_and_strings_do_not_move_scope(self):
        text = 'void sample() { /* } */ const char* s="}"; if (true) { ; } // }\n }'
        lo, hi = S.function_span(text, 'void sample()')
        self.assertEqual(text[lo], '{'); self.assertEqual(hi, len(text)-1)
        with self.assertRaises(ValueError): S.function_span(text, 'void missing()')

    def test_local_finite_identities_do_not_unlock_storage_search(self):
        local_result = {'map_representation': 'finite_physical_descriptor',
                        'finite_error_identity_for_every_event': False,
                        'physical_reference_forcing_retained': True,
                        'all_coefficient_product_graphs_retained': False,
                        'all_configured_branches_bound_to_finite_graph': False}
        with self.assertRaisesRegex(RuntimeError, 'finite-state storage blocked'):
            PLAN.assert_finite_storage_master(local_result)


if __name__ == '__main__': unittest.main()
