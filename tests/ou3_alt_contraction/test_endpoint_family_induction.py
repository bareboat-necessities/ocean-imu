import unittest
from ou3_interval import Interval, matrix_point
from tools.stability.ou3_alt_contraction import endpoint_family_induction as E


class EndpointFamilyInductionTest(unittest.TestCase):
    def test_hull_contains_every_explicit_branch_product(self):
        A1=matrix_point([[0.5,0.0],[0.0,1.0]])
        A2=matrix_point([[0.6,0.0],[0.0,1.0]])
        B1=matrix_point([[0.8,0.0],[0.0,1.0]])
        B2=matrix_point([[0.9,0.0],[0.0,1.0]])
        s0=[E.PartitionSuccessor('p',('root',),E.QualifiedLocalFamily('s0',(A1,A2),True))]
        s1=[E.PartitionSuccessor('q',('p',),E.QualifiedLocalFamily('s1',(B1,B2),True))]
        d=E.prove_endpoint_outer_induction(2,[s0,s1])
        outer=d['endpoint_partitions']['q']
        for B in (B1,B2):
            for A in (A1,A2):
                from ou3_interval import matrix_mul
                self.assertTrue(E.matrix_encloses(outer,matrix_mul(B,A)))
        self.assertTrue(d['outer_union_inclusion_by_induction'])
        self.assertFalse(d['outer_hull_is_realizable_source_history'])

    def test_partitions_prevent_unnecessary_cross_hull(self):
        A=matrix_point([[0.5,0.0],[0.0,1.0]])
        B=matrix_point([[0.9,0.0],[0.0,1.0]])
        step=[
          E.PartitionSuccessor('H',('root',),E.QualifiedLocalFamily('H',(A,),True)),
          E.PartitionSuccessor('A',('root',),E.QualifiedLocalFamily('A',(B,),True)),
        ]
        d=E.prove_endpoint_outer_induction(2,[step])
        self.assertEqual(set(d['endpoint_partitions']),{'H','A'})
        self.assertEqual(d['endpoint_partitions']['H'][0][0].lo,0.5)
        self.assertEqual(d['endpoint_partitions']['A'][0][0].lo,0.9)

    def test_unqualified_or_replay_family_fails_closed(self):
        A=matrix_point([[1.0]])
        with self.assertRaises(ValueError):
            E.advance_partitions(E.initial_partition(1),[
                E.PartitionSuccessor('x',('root',),E.QualifiedLocalFamily('bad',(A,),False))])
        with self.assertRaises(ValueError):
            E.advance_partitions(E.initial_partition(1),[
                E.PartitionSuccessor('x',('root',),E.QualifiedLocalFamily('bad',(A,),True,True))])

    def test_missing_parent_fails_closed(self):
        A=matrix_point([[1.0]])
        with self.assertRaises(ValueError):
            E.advance_partitions(E.initial_partition(1),[
                E.PartitionSuccessor('x',('not-root',),E.QualifiedLocalFamily('x',(A,),True))])

    def test_module_status_requires_actual_600_step_binding(self):
        d=E.build(); self.assertEqual(E.validate(d),[])
        self.assertTrue(d['endpoint_outer_inclusion_induction_materialized'])
        self.assertFalse(d['actual_600_step_source_family_bound'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
