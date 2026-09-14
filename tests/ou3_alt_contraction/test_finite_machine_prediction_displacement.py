"""Full joint24/21-covariance tuner prediction displacement regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_machine_prediction_displacement as X
from tools.stability.ou3_alt_contraction import finite_machine_active_prediction_roots as MPRED
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as EXACT
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
import test_finite_source_bound_live_word as BASE
import test_finite_source_bound_prediction_word as PBASE


def roots(*,changed_sigma=False):
    s=BASE.root_state(); witness,segment,raw,_=PBASE.operands(s)
    physical=SOURCE.QualifiedPhysicalSegment(s.source.root,witness,segment)
    args=PBASE.root_args(); args.pop('temperature_c')
    exact=EXACT.build(s,physical,raw,**args)
    active=s.live.live.live.active
    if changed_sigma:
        sig=tuple(tuple(active.Sigma_aw[i][j]+(F(1,100) if i==j else 0)
                        for j in range(3)) for i in range(3))
        machine_active=ACTIVE.ActiveParameters(active.tau,sig,active.pseudo_period,active.R_S)
    else:
        machine_active=active
    machine=MPRED.build(s,physical,raw,machine_active,mode='separate',
                        exact_active=active,**args)
    return s,segment,raw,exact,machine


class Tests(unittest.TestCase):
    def test_identical_active_parameters_give_exact_zero_full_supply(self):
        s,segment,raw,e,m=roots()
        out=X.compare(s.live.live.live.mekf,segment,raw,e,m,Qbase=s.runtime.Qbase)
        self.assertTrue(all(v==0 for v in out.supply.z))
        self.assertTrue(all(v==0 for row in out.supply.covariance for v in row))
        self.assertEqual(out.exact,out.machine)

    def test_stationary_sigma_change_moves_covariance_not_joint24_mean(self):
        s,segment,raw,e,m=roots(changed_sigma=True)
        out=X.compare(s.live.live.live.mekf,segment,raw,e,m,Qbase=s.runtime.Qbase)
        self.assertTrue(all(v==0 for v in out.supply.z))
        self.assertTrue(any(v!=0 for row in out.supply.covariance for v in row))
        self.assertEqual(out.exact.z,out.machine.z)
        self.assertNotEqual(out.exact.covariance,out.machine.covariance)

    def test_attitude_or_BA_root_discrepancy_cannot_hide_in_tuner_supply(self):
        s,segment,raw,e,m=roots()
        from dataclasses import replace
        # Changing the bias root type/identity must be rejected before paired
        # prediction comparison.  The comparison is reserved for tuner-only
        # coefficient displacement.
        bad=replace(m,bias=replace(m.bias,tau_b=m.bias.tau_b+1))
        with self.assertRaisesRegex(ValueError,'BA-root discrepancy'):
            X.compare(s.live.live.live.mekf,segment,raw,e,bad,Qbase=s.runtime.Qbase)

    def test_readiness_stops_before_postprediction_and_storage(self):
        r=X.readiness()
        self.assertTrue(r['full_prediction_coefficient_displacement_relation_attached'])
        self.assertTrue(r['machine_tau_Sigma_effect_on_full_21x21_covariance_retained'])
        for k in ('post_prediction_floor_scheduler_and_S_service_reexecuted_from_machine_state',
                  'accelerometer_measurement_reexecuted_from_machine_state',
                  'source_uniform_machine_prediction_supply_bound_closed',
                  'machine_root_effect_injected_into_complete_admitted_event_relation',
                  'all_target_libm_and_Eigen_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
