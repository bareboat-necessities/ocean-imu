from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_live_disturbance_overflow as X


class Tests(unittest.TestCase):
    def test_same_quiet_physics_and_finite_source_owned_pulse(self):
        for family in ('BIAS0','BIAS1','BIAS2'):
            before=X.quiet_packet(F(150),bias_family=family)
            pulse=X.quiet_packet(F(150)+X.DT,pulse=True,bias_family=family)
            after=X.quiet_packet(F(150)+2*X.DT,bias_family=family)
            for packet in (before,pulse,after):
                self.assertEqual(packet.physical.acceleration,X.ZERO)
                self.assertEqual(packet.physical.position,X.ZERO)
                self.assertEqual(packet.physical.centered_S,X.ZERO)
                self.assertEqual(packet.physical.gyro_bias,X.ZERO)
                self.assertEqual(packet.internal_accel,X.RAW_ACCEL)
            self.assertEqual(before.internal_gyro,X.ZERO)
            self.assertEqual(after.internal_gyro,X.ZERO)
            self.assertEqual(pulse.internal_gyro,(X.PULSE,0,0))
            with self.assertRaisesRegex(ValueError,'detached from physical rate'):
                replace(pulse,gyro_residual_internal=X.ZERO)

    def test_exact_overflow_and_domain_boundary_remain_explicit(self):
        report=X.build()
        self.assertEqual(X.validate(report),[])
        self.assertGreater(report['squared_component_lower'],report['RNE_finite_overflow_threshold'])
        for key in ('bounded_raw_residual_implies_finite_execution',
                    'nonfinite_execution_has_produced_RestrictedForcing',
                    'target_hardware_counterexample_claimed',
                    'commissioned_startup_contract_falsified',
                    'source_uniform_deployment_supplies_closed','storage_search_allowed'):
            self.assertFalse(report[key])
        report['nonfinite_execution_has_produced_RestrictedForcing']=True
        self.assertTrue(X.validate(report))


class NativeTests(unittest.TestCase):
    def test_ordinary_startup_then_single_finite_pulse_breaks_live_finiteness(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        eigen=Path(os.environ.get('EIGEN_INCLUDE_DIR','/usr/include/eigen3'))
        if not cxx or not (eigen/'Eigen/Dense').is_file():
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1': self.fail('required g++/Eigen unavailable')
            self.skipTest('g++/Eigen unavailable')
        root=Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as td:
            exe=Path(td)/'live-disturbance-overflow'
            subprocess.run([cxx,'-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                '-DEIGEN_DONT_VECTORIZE',f'-I{eigen}',f'-I{root / "src"}',
                str(root/'tests/ou3_alt_contraction/live_disturbance_overflow.cpp'),'-o',str(exe)],
                check=True,capture_output=True,text=True,timeout=180)
            result=subprocess.run([str(exe)],check=True,capture_output=True,text=True,timeout=30)
        rows=[tuple(map(int,line.split())) for line in result.stdout.splitlines()]
        self.assertEqual(rows,[(30002,0,1,1,1,1,1),(30002,1,1,0,0,1,1),
                               (30002,2,1,0,0,1,1),(30002,600,1,0,0,1,1)])


if __name__=='__main__': unittest.main()
