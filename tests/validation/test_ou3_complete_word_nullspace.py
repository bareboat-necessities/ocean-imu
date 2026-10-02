import unittest
from tools.stability.ou3_theorem.complete_word_nullspace import compatibility_line_certificate

class T(unittest.TestCase):
    def test_exact_line(self):
        blocks=[("S",[[1,0,0]]),("mag",[[0,1,-1]])]
        r=compatibility_line_certificate(blocks,[0,1,1])
        self.assertEqual(r["rank"],2); self.assertEqual(r["nullity"],1)
        self.assertTrue(r["kernel_equals_compatibility_line"])

    def test_hidden_second_kernel_fails_closed(self):
        blocks=[("S",[[1,0,0,0]]),("mag",[[0,1,-1,0]])]
        r=compatibility_line_certificate(blocks,[0,1,1,0])
        self.assertEqual(r["nullity"],2)
        self.assertFalse(r["kernel_equals_compatibility_line"])

    def test_wrong_line_fails_closed(self):
        blocks=[("S",[[1,0]]),("mag",[[0,1]])]
        r=compatibility_line_certificate(blocks,[0,1])
        self.assertFalse(r["compatibility_line_annihilated"])
        self.assertFalse(r["kernel_equals_compatibility_line"])

if __name__=="__main__": unittest.main()
