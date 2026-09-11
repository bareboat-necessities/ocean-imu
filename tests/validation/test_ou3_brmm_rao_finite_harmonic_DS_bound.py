import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
MOD=ROOT/'tools'/'stability'/'ou3_brmm_rao_finite_harmonic_DS_bound.py'
spec=importlib.util.spec_from_file_location('ou3_brmm_rao_finite_harmonic_DS_bound',MOD)
M=importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(M)

class RaoFiniteHarmonicDsBoundTests(unittest.TestCase):
    def test_bound_is_derived_and_not_clipped_to_300(self):
        d=M.build()
        self.assertEqual(M.validate(d),[])
        self.assertTrue(d['deterministic_pathwise'])
        self.assertFalse(d['probabilistic_or_PSD_argument_used'])
        self.assertGreater(d['centered_primitive_D_S_upper_m_s'],300.0)
        self.assertFalse(d['legacy_300_m_s_radius_is_proved_sufficient_for_this_bound'])

    def test_subset_does_not_promote_complete_brmm(self):
        d=M.build()
        self.assertFalse(d['complete_BRMM_numeric_qualification_closed'])
        self.assertFalse(d['P4_PASS'])
        self.assertFalse(d['P5_MAY_START'])

if __name__=='__main__':unittest.main()
