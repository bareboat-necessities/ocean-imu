from fractions import Fraction as F
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_startup_timeout_alignment_obstruction as X


class SourceTests(unittest.TestCase):
    def test_docking_source_preserves_circle_without_promoting_infinite_tail(self):
        report=X.docking_source_bounds()
        self.assertLess(report['bounds']['position_norm'],F('8.1'))
        self.assertLess(report['bounds']['direction_mean_if_residual_bounded'],F('.05'))
        self.assertFalse(report['indefinite_tail_residual_bound_proved'])
        self.assertFalse(report['eventual_finite_startup_refuted'])
        for time in (*range(5,150,5),X.DOCKING_TIME):
            yaw,rate=X.docking_yaw_and_rate(time)
            left,left_rate=X.docking_yaw_and_rate(F(time)-F(1,1000000))
            self.assertEqual(left+left_rate/F(1000000),yaw)
            self.assertLessEqual(abs(rate),F('.6'))

    def test_physical_source_and_continuous_direction_budget(self):
        report=X.build()
        self.assertEqual(X.validate(report),[])
        self.assertEqual(F(report['exact']['position_norm_upper']),F('8.09'))
        self.assertEqual(F(report['exact']['velocity_norm_upper']),F('5.1776'))
        self.assertLess(F(report['exact']['direction_mean_norm_upper']),F('.053'))
        self.assertLess(F(report['exact']['direction_primitive_norm_upper']),F('.501'))
        self.assertFalse(report['eventual_finite_startup_refuted'])

    def test_inverted_axis_invariant_is_not_a_reachable_startup_state(self):
        report=X.conditional_inverted_axis_fixed_point()
        self.assertTrue(report['exact_named_binary32_observer_fixed_point'])
        self.assertTrue(report['normalized_projected_accel_z_is_positive_gravity'])
        self.assertEqual(report['quaternion'][0],0)
        self.assertEqual(report['quaternion'][2:],(0,0))
        self.assertFalse(report['source_first_sample_seed_can_be_replaced_by_this_state'])
        self.assertFalse(report['admitted_startup_reaches_inverted_fixed_point'])
        self.assertFalse(report['eventual_finite_startup_refuted'])

    def test_circular_bad_invariant_does_not_claim_docking_or_machine_preservation(self):
        report=X.conditional_circular_antialignment()
        self.assertTrue(report['continuous_observer_and_source_invariant_closed'])
        self.assertGreater(report['projected_force'],0)
        self.assertLess(report['physical_and_temporal_bounds']['position_norm'],F('8.1'))
        self.assertFalse(report['finite_binary32_preservation_qualified'])
        self.assertFalse(report['admitted_startup_reaches_this_invariant'])
        self.assertFalse(report['physical_transition_into_new_circle_qualified'])
        self.assertFalse(report['eventual_finite_startup_refuted'])

    def test_yaw_is_continuous_at_every_rate_change(self):
        for time,expected in ((5,1),(125,F('74.2')),(130,F('76.25')),(135,F('75.25'))):
            value,rate=X.yaw_and_rate(time)
            self.assertEqual(value,expected)
            left,left_rate=X.yaw_and_rate(F(time)-F(1,1000000))
            self.assertEqual(left+left_rate/F(1000000),value)
            self.assertLessEqual(abs(rate),F('.61'))


class NativeTests(unittest.TestCase):
    def test_full_wrapper_reaches_antialigned_docking_and_finite_deadzone_tail(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        eigen=Path(os.environ.get('EIGEN_INCLUDE_DIR','/usr/include/eigen3'))
        if not cxx or not (eigen/'Eigen/Dense').is_file():
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1': self.fail('required g++/Eigen unavailable')
            self.skipTest('g++/Eigen unavailable')
        root=Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as td:
            exe=Path(td)/'docking';packets=Path(td)/'packets.txt'
            subprocess.run([cxx,'-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                '-DEIGEN_DONT_VECTORIZE',f'-I{eigen}',f'-I{root / "src"}',
                str(root/'tests/ou3_alt_contraction/startup_antialigned_docking.cpp'),'-o',str(exe)],
                check=True,capture_output=True,text=True,timeout=180)
            out=subprocess.run([str(exe),str(packets),'40000'],check=True,
                capture_output=True,text=True,timeout=30)
            first,second=map(str.split,out.stdout.splitlines())
            self.assertEqual(first[:2],['dock','30002'])
            self.assertGreater(float(first[3]),2)
            self.assertEqual(first[4:],['live','0'])
            tail=dict(zip(second[::2],second[1::2]))
            self.assertEqual(tail['tail'],'40000')
            for name in ('integral_changed','shadow_differs','live','weight'):
                self.assertEqual(float(tail[name]),0)
            self.assertGreater(float(tail['minlp']),2)
            self.assertLess(float(tail['halferror']),.000001)
            self.assertLess(float(tail['guard']),.03)
            with packets.open() as stream:audit=X.audit_native_packets(stream,docking=True)
            self.assertEqual(audit['samples'],40000)
            self.assertTrue(audit['both_commissioned_profiles_pass'])

    def test_unmodified_wrapper_misses_timeout_and_later_recovers(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        eigen=Path(os.environ.get('EIGEN_INCLUDE_DIR','/usr/include/eigen3'))
        if not cxx or not (eigen/'Eigen/Dense').is_file():
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1': self.fail('required g++/Eigen unavailable')
            self.skipTest('g++/Eigen unavailable')
        root=Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as td:
            exe=Path(td)/'timeout-alignment'; packets=Path(td)/'packets.txt'
            subprocess.run([cxx,'-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                '-DEIGEN_DONT_VECTORIZE',f'-I{eigen}',f'-I{root / "src"}',
                str(root/'tests/ou3_alt_contraction/startup_timeout_alignment_obstruction.cpp'),'-o',str(exe)],
                check=True,capture_output=True,text=True,timeout=180)
            out=subprocess.run([str(exe),str(packets)],check=True,capture_output=True,text=True,timeout=30)
            rows=[line.split() for line in out.stdout.splitlines()]
            self.assertEqual(len(rows),3)
            for row,ordinal in zip(rows[:2],(30002,30602)):
                self.assertEqual(row[:6],['step',str(ordinal),'live','0','proxy','1'])
                self.assertGreater(float(row[7]),1)
                self.assertLess(float(row[9]),.03)
                self.assertEqual(float(row[11]),0)
            self.assertEqual(rows[2][0],'first_live')
            self.assertGreater(int(rows[2][1]),30602)
            self.assertLess(int(rows[2][1]),40000)
            with packets.open() as stream: audit=X.audit_native_packets(stream)
            self.assertEqual(audit['samples'],int(rows[2][1]))
            self.assertTrue(audit['both_commissioned_profiles_pass'])


if __name__=='__main__': unittest.main()
