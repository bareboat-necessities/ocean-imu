from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_startup_disturbance_obstruction as X


class Tests(unittest.TestCase):
    def test_quiet_physics_and_constant_residual_are_one_history(self):
        for family in ('BIAS0','BIAS1','BIAS2'):
            for time in (F(0),F(1,200),F(150),F(10**9)):
                p=X.quiet_packet(time,bias_family=family)
                self.assertEqual(p.physical.acceleration,X.ZERO)
                self.assertEqual(p.physical.position,X.ZERO)
                self.assertEqual(p.physical.centered_S,X.ZERO)
                self.assertEqual(p.physical.beta,X.ZERO)
                self.assertEqual(p.accel_residual_internal,X.RESIDUAL)
                self.assertEqual(p.internal_accel,X.RAW_ACCEL)
                self.assertEqual(tuple(a/X.EPSILON for a in p.internal_accel),(0,0,-1))
        with self.assertRaisesRegex(ValueError,'detached from physical specific force'):
            replace(X.quiet_packet(),accel_residual_internal=X.ZERO)

    def test_literal_guard_and_observer_fixed_points(self):
        first=X.guard_successor(X.GUARD.State())
        second=X.guard_successor(first.state)
        self.assertEqual(replace(second.state,samples=first.state.samples),first.state)
        self.assertEqual(second.conditioned_acc,X.RAW_ACCEL)
        self.assertEqual(second.removed_rms,0)
        self.assertEqual(second.state.weight,0)
        seed=X.observer_fixed_point()
        self.assertEqual(seed.startup.norm,X.EPSILON)
        self.assertEqual(seed.vertical.state,X.VERTICAL.State())
        self.assertLess(seed.startup.norm,X.SEED.rn(F(1,1000)))

    def test_source_margin_is_conditional_and_includes_all_sensor_charges(self):
        args=dict(physical_acceleration_cap=F(44,5),bias_cap=F(2,5),
                  conversion_guard_and_norm_error_cap=F(1,10**5))
        a=X.seed_magnitude_margin(**args,residual_cap=F(3,5))
        b=X.seed_magnitude_margin(**args,residual_cap=F(61,100))
        self.assertGreater(a,0)
        self.assertLess(b,0)
        self.assertEqual(a-b,F(1,100))
        with self.assertRaises(ValueError):
            X.seed_magnitude_margin(**args,residual_cap=-1)

    def test_report_refutes_only_the_arbitrary_bounded_startup_extension(self):
        report=X.build()
        self.assertEqual(X.validate(report),[])
        self.assertTrue(report['finite_prefix_no_initialization_induction'])
        for key in ('arbitrary_bounded_residuals_imply_universal_startup',
                    'post_Live_ISS_bound_supplies_startup_raw_norm_premise',
                    'specified_small_disturbance_startup_theorem_falsified',
                    'hardware_residual_admission_claimed','startup_sensor_contract_selected_by_this_obstruction',
                    'target_libm_compiler_qualification_claimed','storage_search_allowed'):
            self.assertFalse(report[key])
        report['post_Live_ISS_bound_supplies_startup_raw_norm_premise']=True
        self.assertTrue(X.validate(report))


class NativeTests(unittest.TestCase):
    def test_public_wrapper_stays_unseeded_with_or_without_regular_magnetic_service(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        eigen=Path(os.environ.get('EIGEN_INCLUDE_DIR','/usr/include/eigen3'))
        if not cxx or not (eigen/'Eigen/Dense').is_file():
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1': self.fail('required g++/Eigen unavailable')
            self.skipTest('g++/Eigen unavailable')
        root=Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as td:
            exe=Path(td)/'startup-disturbance'
            subprocess.run([cxx,'-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                '-DEIGEN_DONT_VECTORIZE',f'-I{eigen}',f'-I{root / "src"}',
                str(root/'tests/ou3_alt_contraction/startup_disturbance_obstruction.cpp'),'-o',str(exe)],
                check=True,capture_output=True,text=True,timeout=180)
            result=subprocess.run([str(exe)],check=True,capture_output=True,text=True,timeout=30)
        rows=[tuple(map(int,line.split())) for line in result.stdout.splitlines()]
        self.assertEqual(rows,[(0,30602,0,0,0,0,1,1),(1,30602,3825,0,0,1,1,1)])


if __name__=='__main__': unittest.main()
