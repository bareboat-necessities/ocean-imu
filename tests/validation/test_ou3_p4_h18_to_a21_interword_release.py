from __future__ import annotations

import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'/'stability'))

import ou3_p4_h18_to_a21_interword_release as RELEASE


class H18ToA21InterwordReleaseTests(unittest.TestCase):
    def test_shipping_release_topology_and_all_bias_forcing_close_structurally(self):
        d=RELEASE.build()
        self.assertEqual(RELEASE.validate(d),[])
        self.assertTrue(d['release_is_separate_interword_hybrid'])
        self.assertGreaterEqual(d['minimum_complete_H18_words_before_release'],10)
        self.assertTrue(d['all_bias_families_fit_declared_A21_bias_coordinate'])
        self.assertTrue(d['release_structural_source_attachment_closed'])
        self.assertFalse(d['held_bias_error_set_to_zero_at_release'])
        self.assertFalse(d['A21_post_release_basin_landing_closed_here'])
        self.assertFalse(d['P4_PASS'])

    def test_release_maps_are_rectangular_plus_separate_bias_forcing(self):
        J,Eb=RELEASE.release_maps()
        self.assertEqual((len(J),len(J[0])),(21,18))
        self.assertEqual((len(Eb),len(Eb[0])),(21,3))
        for i in range(21):
            for j in range(18):
                expect=1.0 if i==j and i<18 else 0.0
                self.assertEqual((J[i][j].lo,J[i][j].hi),(expect,expect))
        for i in range(21):
            for j in range(3):
                expect=1.0 if i==18+j else 0.0
                self.assertEqual((Eb[i][j].lo,Eb[i][j].hi),(expect,expect))

    def test_every_bias_family_fits_declared_release_coordinate(self):
        d=RELEASE.build()
        self.assertEqual(set(d['all_bias_family_release_bounds']),{'BIAS0','BIAS1','BIAS2'})
        for row in d['all_bias_family_release_bounds'].values():
            self.assertLessEqual(row['held_true_bias_norm_upper_mps2'],row['A21_bias_entry_radius_mps2'])
            self.assertTrue(row['held_bias_forcing_fits_A21_entry_coordinate'])

    def test_false_release_promotion_is_rejected(self):
        d=RELEASE.build();d['held_bias_error_set_to_zero_at_release']=True;d['P4_PASS']=True
        f=RELEASE.validate(d)
        self.assertIn('held_bias_error_set_to_zero_at_release not false',f)
        self.assertIn('P4_PASS not false',f)


if __name__=='__main__':unittest.main()
