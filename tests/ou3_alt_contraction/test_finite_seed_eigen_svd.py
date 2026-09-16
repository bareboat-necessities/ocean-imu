"""Exact scalar graph checks; native comparisons do not promote universality."""
from fractions import Fraction as F
from pathlib import Path
import os
import random
import shutil
import subprocess
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_seed_eigen_svd as X


def inputs():
    b=X.M.round_bits
    cases=[(0,0,b(-1)),(b(F(1,1000)),b(F(2,1000)),b(-1)),
           (b(F(-1,1000)),b(F(-2,1000)),b(-1)),
           (1,2,b(-1)),(b(F(1,100)),0,b(F(-99999,100000))),
           (0,b(F(1,1000)),b(F(-999999,1000000))),
           (b(F(1,1000)),0,b(F(-999999,1000000)))]
    rng=random.Random(528)
    for _ in range(40):
        cases.append((b(F(rng.randrange(-3000,3001),10**6)),
                      b(F(rng.randrange(-3000,3001),10**6)),
                      b(F(-1)+F(rng.randrange(-10,11),10**7))))
    return tuple(cases)


class Tests(unittest.TestCase):
    def test_actual_QR_and_sweeps_compute_the_axis(self):
        pivots=set()
        for v0 in inputs():
            r=X.solve_bits(v0,(0,0,X.ONE)); pivots.add(r.qr.first_pivot)
            self.assertLessEqual(r.sweep_count,3)  # Per-case diagnostic only.
            self.assertEqual(r.axis,tuple(row[2] for row in r.qr.q))
            for v in r.v_history:
                self.assertEqual(r.axis,tuple(row[2] for row in v))
            self.assertFalse(any(X.special(x) for x in r.axis))
        self.assertEqual(pivots,{0,1})

    def test_IEEE_nontrapping_specials_and_zero_signs(self):
        a=X.Arithmetic(); sign=X.M.SIGN
        self.assertEqual(a.div(X.ONE,0),X.INF)
        self.assertEqual(a.div(X.ONE,sign),X.INF|sign)
        self.assertTrue(X.isnan(a.mul(X.INF,0)))
        self.assertTrue(X.isnan(a.sub(X.INF,X.INF)))
        self.assertEqual(a.sqrt(sign),sign)
        self.assertFalse(X.less(X.NAN,X.ONE))
        self.assertFalse(X.less(X.ONE,X.NAN))
        self.assertEqual(a.mul(X.INF,X.M.exact_bits(-1)),X.INF|sign)

    def test_rank_one_axis_and_same_source_witness(self):
        w,r=X.witness((0,0,-1),(0,0,1))
        self.assertEqual(w.axis,(F(1),F(0),F(0)))
        self.assertEqual(w.v0_dot,0)
        self.assertEqual(w.v1_dot,0)
        self.assertEqual(w.norm2_minus_one,0)
        r.validate()
        status=X.readiness()
        self.assertTrue(status['axis_computed_from_scaled_source_vectors_without_solver_input_port'])
        self.assertTrue(status['source_uniform_Jacobi_termination_closed'])
        certificate=status['source_uniform_totality_certificate']
        self.assertEqual(certificate['maximum_Jacobi_sweeps'],2)
        self.assertLess(certificate['second_sweep_residual_upper'],
                        certificate['retained_termination_threshold_lower'])
        self.assertFalse(status['target_compiler_correspondence_closed'])

    def test_source_domain_does_not_cover_arbitrary_solver_vectors(self):
        self.assertTrue(X.in_source_domain((0,0,-1),(0,0,1)))
        self.assertFalse(X.in_source_domain((1,0,-1),(0,0,1)))
        self.assertFalse(X.in_source_domain((0,0,-1),(0,1,0)))
        with self.assertRaisesRegex(ValueError,'source domain'):
            X.witness((1,0,-1),(0,0,1))


HARNESS = r'''
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <cstdint>
#include <cstring>
#include <iostream>
uint32_t bits(float x) { uint32_t b; std::memcpy(&b,&x,4); return b; }
float value(uint32_t b) { float x; std::memcpy(&x,&b,4); return x; }
void emit(float x) { std::cout << ' ' << bits(x); }
int main() {
  uint32_t a,b,c;
  while(std::cin >> a >> b >> c) {
    Eigen::Matrix<float,2,3> m;
    m << value(a),value(b),value(c),0.0f,0.0f,1.0f;
    float scale=m.cwiseAbs().maxCoeff();
    Eigen::Matrix<float,3,2> t=(m/scale).adjoint();
    Eigen::ColPivHouseholderQR<Eigen::Matrix<float,3,2>> qr(t);
    Eigen::Matrix3f v; Eigen::RowVector3f scratch;
    qr.householderQ().evalTo(v,scratch);
    Eigen::Matrix2f work=qr.matrixQR().block(0,0,2,2).triangularView<Eigen::Upper>().adjoint();
    std::cout << qr.colsPermutation().indices()(0);
    emit(scale);
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) emit(v(i,j));
    for(int i=0;i<2;++i) for(int j=0;j<2;++j) emit(work(i,j));
    float maxdiag=work.cwiseAbs().diagonal().maxCoeff();
    unsigned count=0;
    for(;;) {
      const float threshold=std::max(std::numeric_limits<float>::min(),
                                     2.0f*Eigen::NumTraits<float>::epsilon()*maxdiag);
      if(!(std::abs(work(1,0))>threshold || std::abs(work(0,1))>threshold)) break;
      if(++count>32) return 3;
      Eigen::JacobiRotation<float> jl,jr;
      Eigen::internal::real_2x2_jacobi_svd(work,1,0,&jl,&jr);
      work.applyOnTheLeft(1,0,jl); work.applyOnTheRight(1,0,jr);
      v.applyOnTheRight(1,0,jr);
      maxdiag=std::max(maxdiag,std::max(std::abs(work(1,1)),std::abs(work(0,0))));
    }
    std::cout << ' ' << count;
    for(int i=0;i<2;++i) for(int j=0;j<2;++j) emit(work(i,j));
    Eigen::JacobiSVD<Eigen::Matrix<float,2,3>> svd(m,Eigen::ComputeFullV);
    for(int i=0;i<3;++i) emit(svd.matrixV()(i,2));
    std::cout << '\n';
  }
}
'''


class NativeTests(unittest.TestCase):
    def test_pinned_Eigen_QR_work_loop_and_axis_match_bit_for_bit(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        root=Path(os.environ.get('OU3_ALT_TARGET_EIGEN_INCLUDE_DIR',
                  '/tmp/ou3-target/eigen/Eigen-0.3.2/ArduinoEigen'))
        if not cxx or not (root/'Eigen/Dense').is_file():
            self.skipTest('pinned Eigen payload and native C++ compiler required')
        X.REDUCTION.audit_headers(root)
        with tempfile.TemporaryDirectory() as td:
            cpp=Path(td)/'svd.cpp'; exe=Path(td)/'svd'
            cpp.write_text(HARNESS)
            subprocess.run([cxx,'-std=c++20','-O2','-ffp-contract=off','-fno-fast-math',
                            '-DEIGEN_DONT_VECTORIZE',f'-I{root}',str(cpp),'-o',str(exe)],
                           capture_output=True,text=True,check=True,timeout=90)
            rows=subprocess.run([str(exe)],input=''.join(' '.join(map(str,x))+'\n' for x in inputs()),
                                capture_output=True,text=True,check=True,timeout=20).stdout.splitlines()
        self.assertEqual(len(rows),len(inputs()))
        for v0,line in zip(inputs(),rows):
            words=tuple(map(int,line.split())); r=X.solve_bits(v0,(0,0,X.ONE))
            expected=(r.qr.first_pivot,r.qr.scale,
                      *(x for row in r.qr.q for x in row),
                      *(x for row in r.qr.work for x in row),r.sweep_count,
                      *(x for row in r.work_history[-1] for x in row),*r.axis)
            self.assertEqual(words,expected,msg=f'v0={v0}')


if __name__=='__main__': unittest.main()
