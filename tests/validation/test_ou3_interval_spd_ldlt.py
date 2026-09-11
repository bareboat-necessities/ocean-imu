#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
from ou3_interval import Interval,matrix_identity,matrix_mul
import ou3_interval_linear_algebra as L

I=Interval

def contains_identity(A):
    for i,row in enumerate(A):
        for j,x in enumerate(row):
            target=1.0 if i==j else 0.0
            if not x.contains(target):return False
    return True

class TestIntervalSpdLdlt(unittest.TestCase):
    def test_spd_interval_inverse_encloses_identity(self):
        # Symmetric diagonally-dominant interval family with nontrivial widths.
        A=[[I(5.0,5.2),I(-.7,.7),I(-.4,.4)],
           [I(-.7,.7),I(4.5,4.8),I(-.6,.6)],
           [I(-.4,.4),I(-.6,.6),I(3.8,4.1)]]
        X=L.matrix_inverse_spd_ldlt(A)
        self.assertTrue(contains_identity(matrix_mul(A,X)))
        Y=L.matrix_inverse_gauss_jordan(A)
        self.assertTrue(contains_identity(matrix_mul(A,Y)))

    def test_indefinite_symmetric_family_not_accepted_by_ldlt(self):
        A=[[I(-1,1),I(0,0)],[I(0,0),I(1,2)]]
        with self.assertRaises(L.IntervalPivotError):
            L.matrix_inverse_spd_ldlt(A)

    def test_nonsymmetric_general_inverse_path_still_available(self):
        A=[[I(2,2),I(1,1)],[I(0,0),I(3,3)]]
        X=L.matrix_inverse_gauss_jordan(A)
        self.assertTrue(contains_identity(matrix_mul(A,X)))

if __name__=='__main__':unittest.main()
