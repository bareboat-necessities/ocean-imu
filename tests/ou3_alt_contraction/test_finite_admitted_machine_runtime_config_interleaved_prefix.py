"""Machine guard/tracker-LPF runtime-configuration ancestry regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_runtime_config_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
import test_finite_admitted_machine_vertical_stillness_interleaved_prefix as BASE
import test_finite_source_bound_live_word as LBASE


class Tests(unittest.TestCase):
    def test_default_fixture_binds_guard_and_tracker_cutoff_to_shipping_runtime(self):
        lower,_,_,_,_,_,_=BASE.executed_fixture()
        s=X.begin(lower)
        self.assertEqual(s.base.guard_cfg,X._machine_guard_cfg(X._runtime(s.base)))
        self.assertEqual(s.base.separate_source.lpf.cutoff_hz,X.DEFAULT_TRACKER_CUTOFF)
        self.assertEqual(s.base.fma_source.lpf.cutoff_hz,X.DEFAULT_TRACKER_CUTOFF)
        self.assertEqual(s.config_steps,0)

    def test_detached_machine_guard_config_is_rejected_at_entry(self):
        lower,_,_,_,_,_,_=BASE.executed_fixture()
        bad=replace(lower,guard_cfg=replace(lower.guard_cfg,cutoff_hz=B.rn32(13)))
        with self.assertRaisesRegex(ValueError,'guard configuration detached'):
            X.begin(bad)

    def test_tracker_cutoff_cannot_be_independently_spliced(self):
        lower,_,_,_,_,_,_=BASE.executed_fixture()
        lp=replace(lower.separate_source.lpf,cutoff_hz=B.rn32(5))
        src=replace(lower.separate_source,lpf=lp)
        bad=replace(lower,separate_source=src)
        with self.assertRaisesRegex(ValueError,'tracker LPF cutoff detached'):
            X.begin(bad)

    def test_same_admitted_event_advances_runtime_config_validation_count(self):
        lower,kwargs,machine,witness,_,_,_=BASE.executed_fixture()
        s=X.begin(lower)
        out=X.imu_step(s,**witness,**machine,**kwargs)
        self.assertEqual(out.state.config_steps,1)
        self.assertEqual(out.state.base.source_steps,s.base.source_steps+1)
        self.assertEqual(out.state.base.guard_cfg,X._machine_guard_cfg(X._runtime(out.state.base)))

    def test_MAG_and_HOLD_preserve_configuration_validation_ordinal(self):
        lower,_,_,_,_,_,_=BASE.executed_fixture(); s=X.begin(lower); n=s.config_steps
        word=lower.base.base.prefix.prefix.live.live_word
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word)); self.assertEqual(s.config_steps,n)
        s,_=X.set_hold(s,hold=False); self.assertEqual(s.config_steps,n)

    def test_complete_word_cannot_skip_configuration_checks(self):
        lower,_,_,_,_,_,_=BASE.executed_fixture()
        with self.assertRaises((ValueError,TypeError)):
            X.complete(X.begin(lower))

    def test_readiness_closes_current_word_config_not_mutable_setters_or_storage(self):
        r=X.readiness()
        for k in ('machine_guard_config_binary32_projection_bound_to_carried_RuntimeConfig',
                  'machine_guard_runtime_config_ancestry_closed_for_current_word',
                  'tracker_LPF_shipping_default_constructor_reset_source_shape_checked',
                  'tracker_LPF_cutoff_bound_to_default_constructor_reset_for_current_word',
                  'separate_and_FMA_tracker_LPF_cutoff_cannot_diverge',
                  'complete_word_requires_runtime_config_check_on_all_600_IMU_edges'):
            self.assertTrue(r[k])
        for k in ('tracker_LPF_mutable_setter_ancestry_closed',
                  'private_Mahony_mutable_config_setter_ancestry_closed',
                  'startup_machine_runtime_config_history_attached',
                  'target_libm_Eigen_and_compiler_profile_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
