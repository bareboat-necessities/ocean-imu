from __future__ import annotations
import math,sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.rank_loss_interval_factor import (
    FactorState,exact,eye,zeros,verified_inverse,source_covariance,
    schur_scalar_information,residualized_gram,generalized_ratio_lower,
    source_range_audit,
)

class RankLossIntervalFactorTests(unittest.TestCase):
    def test_literal_four_s_point_box_certifies(self):
        from tools.stability.ou3_theorem.rank_loss_literal_boxes import four_s_gamma_box
        cells=((.10,.10),(.35,.35),(.65,.65),(.95,.95))
        z=four_s_gamma_box(cells,(1.0,1.0))
        self.assertTrue(z["verified"])
        self.assertGreater(z["lower"],0.0)

    def test_literal_source_certificate_fails_closed_on_magnetic(self):
        from tools.stability.ou3_theorem.rank_loss_literal_boxes import rank_loss_factor_certificate
        z=rank_loss_factor_certificate(max_depth=2)
        self.assertFalse(z["source_uniform_verified"])
        self.assertEqual(z["gamma_M_lower"],0.0)
        self.assertEqual(z["beta_lower"],0.0)

    def test_verified_generic_inverse(self):
        a=exact([[2.0,.1],[.1,1.0]])
        ai,c=verified_inverse(a)
        self.assertTrue(c["verified"])
        self.assertLess(c["residual_ratio_upper"],1e-12)
        self.assertAlmostEqual(ai.mid[0][0],1/1.99,places=12)

    def test_factor_propagation_preserves_shared_source(self):
        st=FactorState(eye(2),zeros(2,0))
        st=st.predict(exact([[1,1],[0,1]]),exact([[1],[0]]))
        o,c=st.observe(exact([[1,0]]),exact([[2]]))
        self.assertEqual(o.shape,(1,2));self.assertEqual(c.shape,(1,2))
        self.assertAlmostEqual(c.mid[0][0],1.0)
        self.assertAlmostEqual(c.mid[0][1],2.0)

    def test_four_row_schur_scalar_exact_case(self):
        # va is orthogonal to span(V0), R=I => gamma=||va||^2=1.
        v0=exact([[1,0,0],[0,1,0],[0,0,1],[0,0,0]])
        va=exact([[0],[0],[0],[1]])
        z=schur_scalar_information(v0,va,eye(4))
        self.assertTrue(z["verified"]);self.assertGreater(z["lower"],.999999999)

    def test_residualized_magnetic_gram(self):
        obs=exact([[1,0],[0,1],[0,0]])
        nuis=exact([[0],[0],[1]])
        g,z=residualized_gram(obs,nuis,eye(3))
        self.assertTrue(z["verified"]);self.assertGreater(z["lower"],.999999999)
        self.assertEqual(g.shape,(2,2))

    def test_generalized_ratio_lower(self):
        z=generalized_ratio_lower(exact([[2,0],[0,3]]),exact([[4,0],[0,5]]))
        self.assertTrue(z["verified"]);self.assertGreaterEqual(z["lower"],.4-1e-15)

    def test_source_range_audit_fails_closed(self):
        z=source_range_audit()
        self.assertFalse(z["source_uniform_verified"])
        self.assertEqual(z["gamma_S_lower"],0.0)
        self.assertEqual(z["gamma_M_lower"],0.0)
        self.assertEqual(z["beta_lower"],0.0)

if __name__=="__main__":unittest.main()
