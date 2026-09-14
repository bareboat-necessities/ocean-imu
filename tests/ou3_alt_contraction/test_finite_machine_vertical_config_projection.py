"""Regression for the machine private-Mahony runtime-config projection."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_vertical_stillness_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
import test_finite_admitted_machine_tunestate_interleaved_prefix as TBASE


class Tests(unittest.TestCase):
    def test_projection_uses_settle_sec_not_an_unrelated_guard_field(self):
        mt=TBASE.state(usable=True)
        runtime=X._runtime(mt)
        got=X._machine_vertical_cfg(runtime)
        src=runtime.vertical_cfg
        self.assertEqual(got.two_kp,B.rn32(src.two_kp))
        self.assertEqual(got.two_ki,B.rn32(src.two_ki))
        self.assertEqual(got.gravity,B.rn32(src.gravity))
        self.assertEqual(got.settle_sec,B.rn32(src.settle_sec))

    def test_settle_sec_change_changes_only_projected_settle_sec(self):
        mt=TBASE.state(usable=True)
        runtime=X._runtime(mt); src=runtime.vertical_cfg
        changed=replace(runtime,vertical_cfg=replace(src,settle_sec=F(src.settle_sec)+F(1,8)))
        a=X._machine_vertical_cfg(runtime); b=X._machine_vertical_cfg(changed)
        self.assertEqual((a.two_kp,a.two_ki,a.gravity),(b.two_kp,b.two_ki,b.gravity))
        self.assertNotEqual(a.settle_sec,b.settle_sec)
        self.assertEqual(b.settle_sec,B.rn32(changed.vertical_cfg.settle_sec))


if __name__=='__main__': unittest.main()
