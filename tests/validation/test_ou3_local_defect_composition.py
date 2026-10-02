import unittest
from tools.stability.ou3_theorem.local_defect_composition import certificate, verify_endpoint, event_boundary_pairing, left_attitude_error, carried_error, literal_mean_factor

class T(unittest.TestCase):
    def test_exact_local_composition(self):
        r=certificate()
        self.assertTrue(r["synthetic_endpoint_identity_exact"])
        self.assertFalse(r["endpoint_residual_used_as_input"])
        self.assertFalse(r["native_literal_boundary_export_complete"])

    def test_left_attitude_error_sign(self):
        import math
        a=.02
        # estimated W->B identity; true B->W=Rx(+a), hence true W->B=Rx(-a)
        true_bw=[math.sin(a/2),0,0,math.cos(a/2)]
        e=left_attitude_error([0,0,0,1],true_bw)
        self.assertAlmostEqual(e[0],-a,places=12)
        self.assertAlmostEqual(e[1],0,places=12)
        self.assertAlmostEqual(e[2],0,places=12)

    def test_carried_error_layout(self):
        e={"estimator_state":[[0]]*21,"estimator_quaternion":[[0],[0],[0],[1]],
           "physical_quaternion":[[0],[0],[0],[1]],"physical_bg":[[0],[0],[0]],
           "physical_v":[[1],[2],[3]],"physical_p":[[4],[5],[6]],
           "physical_S":[[7],[8],[9]],"physical_a":[[10],[11],[12]],
           "physical_ba":[[13],[14],[15]]}
        x=carried_error(e)
        self.assertEqual(len(x),21)
        self.assertEqual([float(x[i][0]) for i in range(6,21)],list(range(1,16)))

    def test_boundary_pairing_skips_only_mean_neutral_events(self):
        ev=[{"kind":"prediction"},{"kind":"sync"},{"kind":"sync_completion"},
            {"kind":"correction"},{"kind":"reset"}]
        self.assertEqual(event_boundary_pairing(ev),[(0,3),(3,4)])

    def test_nontrivial_three_step(self):
        b=[
          {"A":[["1"]],"e_before":[["2"]],"e_after":[["5"]]},
          {"A":[["2"]],"e_before":[["5"]],"e_after":[["9"]]},
          {"A":[["3"]],"e_before":[["9"]],"e_after":[["28"]]}]
        self.assertTrue(verify_endpoint(b)["endpoint_identity_exact"])

if __name__=="__main__":
    unittest.main()
