"""Target-operation branch inclusion and persistent-history attachment."""
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
import unittest

from tools.stability.ou3_alt_contraction import target_wpe_compiler as T
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as R


def q(x): return T.B.rn32(F(x))


def oracle(kind,arg):
    # A checked test witness, not an assertion about host/target libm accuracy.
    with localcontext() as ctx:
        ctx.prec=80
        x=Decimal(arg.numerator)/Decimal(arg.denominator)
        result={'exp':x.exp,'log':x.ln,'sqrt':x.sqrt}[kind]()
        return T.B.rn32(F(result))


class Tests(unittest.TestCase):
    def test_reset_target_step_lifts_into_same_FMA_product(self):
        before=T.initial()
        selected=T.step(before,dt=T.U.DT,vertical_accel=q(F(1,8)),libm=oracle)
        self.assertEqual(selected.moment.branch,'pre-moment-start')
        cfg=R.WPEConfig(T.U.LAMBDA,4,F(1,20),20,180)
        product=T.X.initial(cfg,libm_profile=T.U.ERROR_PROFILE)
        after,_=T.attach_to_product(product,selected,separate=selected.witnesses)
        self.assertEqual(after.fma,selected.state.moment)
        self.assertEqual(after.logs.fma,selected.state.log)
        with self.assertRaisesRegex(ValueError,'pre-state'):
            T.attach_to_product(replace(product,fma=replace(product.fma,hp1=q(1))),
                                selected,separate=selected.witnesses)

    def test_each_post_start_early_branch_has_a_selected_target_successor(self):
        scenarios=(
            (T.M.State(elapsed=q(40)),'insufficient-weight'),
            (T.M.State(elapsed=q(40),weight=q(1)),'degenerate-moments'),
            (T.M.State(elapsed=q(40),weight=q(1),velocity_sq=q(F(1,1000)),
                       elevation_sq=q(1)),'nonpositive-omega'),
        )
        for moment,expected in scenarios:
            with self.subTest(expected=expected):
                result=T.step(T.State(moment),dt=T.U.DT,vertical_accel=q(0),libm=oracle)
                self.assertEqual(result.moment.branch,expected)
                self.assertIsNone(result.log_step)
                self.assertIsNone(result.witnesses.raw_log)

    def test_raw_period_and_fused_log_update_share_actual_target_moments(self):
        moment=T.M.State(elapsed=q(40),weight=q(1),velocity_sq=q(9),elevation_sq=q(1))
        initial=T.step(T.State(moment),dt=T.U.DT,vertical_accel=q(0),libm=oracle)
        self.assertEqual(initial.moment.branch,'valid-period')
        self.assertEqual(initial.witnesses.raw_log.raw_period,initial.moment.raw_period)
        self.assertTrue(initial.log_step.contracted)
        next_step=T.step(initial.state,dt=T.U.DT,vertical_accel=q(F(1,4)),libm=oracle)
        self.assertEqual(next_step.moment.branch,'valid-period')
        self.assertTrue(next_step.log_step.contracted)
        self.assertEqual(next_step.witnesses.log.sea_period_exp,
                         next_step.witnesses.horizon.period_exp)
        self.assertEqual(next_step.log_step.before,initial.state.log)

    def test_compiler_mapping_rejects_removed_fused_site(self):
        parts=['.file 2 "/ocean-imu/src/tuner/WavePeriodEstimator.h"']
        for line,ops in T.EXPECTED.items():
            parts.append('.loc 2 '+str(line))
            parts.extend(op+' f0, f1, f2' for op in ops)
        parts.append('\t.size\t_ZN19WavePeriodEstimator6updateEff')
        valid='\n'.join(parts)
        T.arithmetic_sites(valid)
        with self.assertRaisesRegex(ValueError,'operation selection changed'):
            T.arithmetic_sites(valid.replace('msub.s','sub.s',1))

    def test_named_target_profile_preserves_whole_firmware_obligation(self):
        r=T.readiness()
        self.assertTrue(r['pinned_WPE_compiler_profile_selection_closed'])
        self.assertEqual(r['persistent_actual_target_track'],'fma')
        self.assertFalse(r['whole_firmware_compiler_and_link_correspondence_closed'])

    def test_selection_rejects_a_different_compiler_or_flags(self):
        path=T.Path(T.__file__).resolve().parents[1]/'ou3_alt_target_wpe_compiler.json'
        report=json.loads(path.read_text())
        for key,changed in (('compiler_driver_sha256','0'*64),
                            ('sdk_cpp_flags_sha256','1'*64),
                            ('compiler_flags',report['compiler_flags']+['-ffast-math'])):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError,'pinned target build'):
                T.readiness(report={**report,key:changed})


if __name__=='__main__': unittest.main()
