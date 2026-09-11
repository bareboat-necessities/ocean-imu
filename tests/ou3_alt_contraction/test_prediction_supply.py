import unittest
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import prediction_supply as P

class PredictionSupplyTests(unittest.TestCase):
    def test_exact_block_metric_supply(self):
        d=P.build();self.assertEqual(P.validate(d),[])
        self.assertTrue(d['prediction_supply_structure_closed'])
        self.assertTrue(d['exact_rational_regression']['inequality_exact'])
        self.assertFalse(d['F_invertibility_required'])
        self.assertFalse(d['Young_factor_or_epsilon_used'])

    def test_shipping_process_Q_is_strict_in_both_modes(self):
        d=P.build();self.assertGreater(d['H18_shipping_Q_lambda_min_lower'],0);self.assertGreater(d['A21_shipping_Q_lambda_min_lower'],0)
        self.assertFalse(d['source_uniform_prediction_forcing_bound_closed_here'])

if __name__=='__main__':unittest.main()
