import unittest
from tools.stability.ou3_theorem.accel_conditional_information import certificate
class AccelConditionalInformationTests(unittest.TestCase):
 def test_schur_identity(self):
  c=certificate();self.assertLess(c["schur_identity_defect"],1e-10);self.assertFalse(c["shipping_service_proved"])
if __name__=="__main__":unittest.main()
