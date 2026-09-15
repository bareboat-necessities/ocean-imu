import unittest
from tools.stability.ou3_alt_contraction import admitted_startup_mahony_invariant as M

class Tests(unittest.TestCase):
    def test_shared_seed_and_binary32_invariant_close_for_alt(self):
        x=M.build()
        self.assertEqual(M.validate(x),[])
        self.assertTrue(x['admitted_source_private_Mahony_startup_to_Live_invariant_closed'])
        self.assertFalse(x['magnetic_accumulation_frame_accuracy_closed_by_this_invariant'])
        self.assertFalse(x['storage_search_allowed'])
        self.assertFalse(x['commissioned_startup_sensor_profiles_covered'])
        x['commissioned_startup_sensor_profiles_covered']=True
        self.assertIn('commissioned_startup_sensor_profiles_covered not false',M.validate(x))

if __name__=='__main__': unittest.main()
