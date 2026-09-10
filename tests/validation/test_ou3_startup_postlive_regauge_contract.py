import unittest

from tools.stability import ou3_startup_postlive_regauge_contract as C


class StartupPostLiveRegaugeContractTest(unittest.TestCase):
    def test_contract(self):
        d = C.build()
        self.assertEqual([], C.validate(d))
        self.assertTrue(d["STRUCTURAL_POSTLIVE_REGAUGE_PATH_CLOSED"])
        self.assertFalse(d["startup_magnetic_PE_required_for_timeout_branch"])
        self.assertFalse(d["P5_PASS"])


if __name__ == "__main__":
    unittest.main()
