from fractions import Fraction as F
import unittest
from tools.stability.ou3_theorem.paired_two_abel import *
class T(unittest.TestCase):
 def test_boundary(self): self.assertEqual(collapse_slab_boundaries([([1,0],[-2,3]),([2,-3],[-4,5])]),[[F(1),F(0)],[F(0),F(0)],[F(-4),F(5)]])
 def test_pair(self): self.assertEqual(frobenius_sq(paired_source_columns([[1,2]],[[-1,-2]])),0)
 def test_open(self): self.assertIsNone(certificate_schema()["source_uniform_induced_norm_numeric"])
if __name__=="__main__": unittest.main()
