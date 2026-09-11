import dataclasses
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS


class SourceUniformBiasPrefixLineageTests(unittest.TestCase):
    def _lineage(self, n=3):
        out=[]
        parent="root"
        for i in range(n):
            child=f"c{i+1}"
            out.append(SELECTORS.PrefixSelector(
                sample_index=i,prefix_length=i+1,parent_branch_ordinal=0,successor_ordinal=0,
                parent_source_cell_id=parent,source_cell_id=child,
                sample_coordinates=None,active_schedule=None,actual_rs_std_xyz=None,
                H_before=None,H_after=None,A_before=None,A_after=None,
                H_events_this_sample=(),A_events_this_sample=(),H_event_cells=(),A_event_cells=(),
                H_floor_case=None,A_floor_case=None,
            ))
            parent=child
        return out

    def test_certificate_validates_without_promoting_p4(self):
        d=BIAS.build()
        self.assertEqual(BIAS.validate(d),[])
        self.assertTrue(d['source_uniform_projection_bias_coordinate_materialized'])
        self.assertTrue(d['joint_shared_w_supply_still_required'])
        self.assertFalse(d['independent_per_sample_absolute_bias_boxes_used'])
        self.assertFalse(d['production_selector_lineages_attached_here'])
        self.assertFalse(d['P4_PASS'])
        self.assertEqual(d['P3_delta'],1e-18)

    def test_every_family_propagates_one_recursive_prefix_state(self):
        lineage=self._lineage(4)
        for family in BIAS.FAMILIES:
            x=BIAS.attach_to_selector_lineage(family,lineage)
            self.assertTrue(x.source_uniform)
            self.assertTrue(x.recurrence_outer_relation)
            self.assertEqual(x.endpoint_source_cell_id,'c4')
            self.assertEqual(set(x.prefix_boxes),{'c1','c2','c3','c4'})
            for cell in x.prefix_boxes.values():
                for axis in cell:
                    self.assertGreaterEqual(axis.lo,-x.absolute_component_cap)
                    self.assertLessEqual(axis.hi,x.absolute_component_cap)

    def test_driver_or_factor_mutation_changes_recursive_successor(self):
        p=BIAS.family_parameters('BIAS1')
        cap=p['cap']
        root=BIAS.Interval(-cap/4,cap/2)
        nominal=BIAS._intersect(p['phi']*root+p['driver'],BIAS.Interval(-cap,cap))
        no_driver=BIAS._intersect(p['phi']*root,BIAS.Interval(-cap,cap))
        self.assertNotEqual((nominal.lo,nominal.hi),(no_driver.lo,no_driver.hi))

    def test_broken_selector_ancestry_rejected(self):
        lineage=self._lineage(3)
        lineage[2]=dataclasses.replace(lineage[2],parent_source_cell_id='detached')
        with self.assertRaisesRegex(ValueError,'ancestry broken'):
            BIAS.attach_to_selector_lineage('BIAS0',lineage)

    def test_families_are_not_collapsed(self):
        d=BIAS.build()
        self.assertEqual(set(d['family_prefix_relations']),set(BIAS.FAMILIES))
        self.assertFalse(d['three_bias_families_collapsed_into_generic_box'])


if __name__=='__main__':
    unittest.main()
