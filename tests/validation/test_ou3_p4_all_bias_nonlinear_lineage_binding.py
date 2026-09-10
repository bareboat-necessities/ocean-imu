import dataclasses
import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'/'stability'))

from ou3_interval import Interval
import ou3_p4_all_bias_nonlinear_lineage_binding as BIND
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS

I=Interval.point


class DummyConstants:
    accel_bias_tau_s=I(3600.0)


def lineage(n=3):
    out=[]; parent='root'
    for k in range(n):
        child=f'c{k+1}'
        out.append(SELECTORS.PrefixSelector(
            sample_index=k,prefix_length=k+1,parent_branch_ordinal=0,successor_ordinal=0,
            parent_source_cell_id=parent,source_cell_id=child,sample_coordinates=None,
            active_schedule=None,actual_rs_std_xyz=None,H_before=None,H_after=None,
            A_before=None,A_after=None,H_events_this_sample=(),A_events_this_sample=(),
            H_event_cells=(),A_event_cells=(),H_floor_case=None,A_floor_case=None))
        parent=child
    return out


class AllBiasNonlinearLineageBindingTests(unittest.TestCase):
    def test_status_is_all_family_and_fail_closed(self):
        d=BIND.build()
        self.assertEqual(BIND.validate(d),[])
        self.assertEqual(d['required_bias_families'],['BIAS0','BIAS1','BIAS2'])
        self.assertTrue(d['all_three_bias_families_have_executable_A21_lineage_binding'])
        self.assertTrue(d['joint_shared_w_tau_mismatch_supply_retained_for_augmented_master'])
        self.assertFalse(d['homogeneous_BIAS1_assumption_used_for_all_families'])
        self.assertFalse(d['P4_PASS'])

    def test_each_family_adapter_accepts_its_own_recursive_prefix_certificate(self):
        xs=lineage()
        for family in ('BIAS0','BIAS1','BIAS2'):
            a=BIND.adapter_for_lineage(family,xs)
            self.assertTrue(a.source_uniform_materialization)
            a.validate_homogeneous(xs,DummyConstants(),0.4)
            self.assertEqual(a.endpoint_source_cell_id,'c3')

    def test_detached_prefix_bias_cell_is_rejected(self):
        xs=lineage()
        a=BIND.adapter_for_lineage('BIAS2',xs)
        bad=dict(a.bias_true_by_source_cell_id)
        bad['c2']=(I(99.0),I(99.0),I(99.0))
        mutated=dataclasses.replace(a,bias_true_by_source_cell_id=bad)
        with self.assertRaisesRegex(RuntimeError,'detached from recursive prefix certificate'):
            mutated.validate_homogeneous(xs,DummyConstants(),0.4)

    def test_missing_prefix_is_rejected(self):
        xs=lineage()
        a=BIND.adapter_for_lineage('BIAS0',xs)
        bad=dict(a.bias_true_by_source_cell_id);bad.pop('c2')
        mutated=dataclasses.replace(a,bias_true_by_source_cell_id=bad)
        with self.assertRaisesRegex(RuntimeError,'every and only selector prefix'):
            mutated.validate_homogeneous(xs,DummyConstants(),0.4)


if __name__=='__main__':unittest.main()
