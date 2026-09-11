"""Exact covariance-prediction identities; not source admission or stability evidence."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_covariance as C

def sym_cov(n=21):
    u=[F((i%5)-2,17) for i in range(n)]
    return M.plus(M.scaled(M.eye(n),F(3,2)),[[x*y for y in u] for x in u])

def blocks(active=True,corr=True):
    faa=M.eye(6); faa[0][3]=faa[1][4]=faa[2][5]=F(1,200); qaa=M.scaled(M.eye(6),F(1,100000))
    coeff=[(F(1,201),F(1,81000),F(1,48100000),F(99,100),F(1,200))]*3; fll=C.linear_transition(coeff)
    if corr:
        qu=M.scaled(M.eye(4),F(1,10000)); qu[0][3]=qu[3][0]=F(1,20000); sig=M.eye(3); sig[0][1]=sig[1][0]=F(1,10)
        qll=C.correlated_linear_process(qu,sig)
    else: qll=C.independent_linear_process([M.scaled(M.eye(4),F(i+1,10000)) for i in range(3)])
    ph=F(199,200) if active else F(1); qbb=M.scaled(M.eye(3),F(1,1000000)) if active else M.zeros(3,3)
    return C.Blocks(faa,qaa,fll,qll,ph,qbb,active)

class Tests(unittest.TestCase):
    def test_shipping_block_order_equals_dense_congruence(self):
        for active in (False,True):
            for corr in (False,True):
                with self.subTest(active=active,corr=corr):
                    b=blocks(active,corr); p=sym_cov(); self.assertEqual(C.shipping_block_successor(p,b),C.dense_successor(p,b)); self.assertEqual(C.validate_shipping_identity(p,b),C.dense_successor(p,b))
    def test_linear_transition_is_same_mean_coefficient_map(self):
        coeff=[(F(2,5),F(3,7),F(4,9),F(5,11),F(1,13))]*3; A=C.linear_transition(coeff); x=[F(i-4,10) for i in range(12)]; y=M.mv(A,x)
        for a in range(3):
            v,p,S,w=x[a],x[3+a],x[6+a],x[9+a]; va,pa,Sa,alpha,h=coeff[a]
            self.assertEqual((y[a],y[3+a],y[6+a],y[9+a]),(v+va*w,p+h*v+pa*w,S+h*p+h*h*v/2+Sa*w,alpha*w))
    def test_correlated_Q_uses_same_sigma_in_every_group_block(self):
        q=[[F(i+j+1,19) for j in range(4)] for i in range(4)]; q=[[q[i][j] if i<=j else q[j][i] for j in range(4)] for i in range(4)]; sig=[[F(2,3),F(1,7),0],[F(1,7),F(3,4),F(1,11)],[0,F(1,11),F(4,5)]]; out=C.correlated_linear_process(q,sig); off=(0,3,6,9)
        for g in range(4):
            for h in range(4):
                for i in range(3):
                    for j in range(3): self.assertEqual(out[off[g]+i][off[h]+j],q[g][h]*sig[i][j])
    def test_independent_Q_has_no_cross_axis_terms(self):
        out=C.independent_linear_process([M.scaled(M.eye(4),F(k+1,13)) for k in range(3)]); idx=(0,3,6,9)
        for a in range(3):
            for b in range(3):
                if a!=b:
                    for i in range(4):
                        for j in range(4): self.assertEqual(out[idx[i]+a][idx[j]+b],0)
    def test_held_BA_requires_literal_identity_branch(self):
        b=blocks(False); self.assertEqual(b.phi_b,1); self.assertEqual(list(map(list,b.Q_BB)),M.zeros(3,3))
        with self.assertRaises(ValueError): C.Blocks(b.F_AA,b.Q_AA,b.F_LL,b.Q_LL,F(99,100),b.Q_BB,False)
        with self.assertRaises(ValueError): C.Blocks(b.F_AA,b.Q_AA,b.F_LL,b.Q_LL,1,M.eye(3),False)
    def test_active_BA_cross_blocks_share_exact_same_phi(self):
        b=blocks(True); p=sym_cov(); out=C.shipping_block_successor(p,b); FA=list(map(list,b.F_AA)); FL=list(map(list,b.F_LL)); ph=b.phi_b
        self.assertEqual([r[18:] for r in out[:6]],M.scaled(M.mm(FA,[r[18:] for r in p[:6]]),ph)); self.assertEqual([r[18:] for r in out[6:18]],M.scaled(M.mm(FL,[r[18:] for r in p[6:18]]),ph))
    def test_rejects_asymmetric_process_blocks(self):
        b=blocks(True); bad=[list(r) for r in b.Q_LL]; bad[0][1]+=1
        with self.assertRaises(ValueError): C.Blocks(b.F_AA,b.Q_AA,b.F_LL,bad,b.phi_b,b.Q_BB,True)
    def test_readiness_is_fail_closed(self):
        r=C.readiness(); self.assertTrue(r['shipping_covariance_block_identity']); self.assertFalse(r['attitude_F_Q_primitives_source_attached']); self.assertFalse(r['pending_aw_covariance_inflation_attached']); self.assertFalse(r['periodic_S_service_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
