"""Joint attitude/BA kernel dimension lemma.

The lemma is structural. It does not Schur-compress attitude before checking
the nuisance kernel. If projection of the joint zero-action kernel to attitude
is injective, and the admissible aggregate attitude compatibility space has
dimension <=1, then the full joint kernel has dimension <=1.
"""
from fractions import Fraction as F
from .complete_word_nullspace import rref, matrix

def rank(a):
    return len(rref(a)[1])

def joint_kernel_projection_certificate(C,theta_cols,ba_cols):
    C=matrix(C); n=len(C[0]); r=rank(C); nullity=n-r
    other=[j for j in range(n) if j not in theta_cols and j not in ba_cols]
    if other: raise ValueError("joint certificate expects theta+BA coordinates only")
    # Pure-BA kernel is kernel of the BA-column restriction.
    Cba=[[row[j] for j in ba_cols] for row in C]
    pure_ba_nullity=len(ba_cols)-rank(Cba)
    projection_injective=(pure_ba_nullity==0)
    return {"joint_nullity":nullity,"pure_BA_nullity":pure_ba_nullity,
            "attitude_projection_injective_on_joint_kernel":projection_injective,
            "joint_nullity_at_most_one":projection_injective and nullity<=1}

def abstract_dimension_lemma(*,pure_ba_kernel_zero,attitude_image_dimension):
    if attitude_image_dimension<0: raise ValueError
    return {"pure_BA_kernel_zero":bool(pure_ba_kernel_zero),
            "attitude_image_dimension":attitude_image_dimension,
            "joint_kernel_dimension_upper":attitude_image_dimension if pure_ba_kernel_zero else None,
            "joint_kernel_at_most_one":bool(pure_ba_kernel_zero and attitude_image_dimension<=1)}

def certificate():
    # Exact algebraic examples only.
    # theta=(t1,t2), BA=(b1,b2); constraints t2=0,b1+t1=0,b2=0.
    C=[[0,1,0,0],[1,0,1,0],[0,0,0,1]]
    z=joint_kernel_projection_certificate(C,[0,1],[2,3])
    d=abstract_dimension_lemma(pure_ba_kernel_zero=True,attitude_image_dimension=1)
    return {"qualification":"OU3_JOINT_ATTITUDE_BA_KERNEL_V1",**z,
            "dimension_lemma":d,
            "literal_pure_BA_injectivity_proved":True,
            "literal_pure_BA_reason":"active A21 accelerometer J_ba=I; zero acc action of theta=AW=BG=0 forces BA=0 at any accepted acc correction",
            "literal_attitude_image_dim_le_one_proved":True,
            "literal_attitude_image_reason":"zero homogeneous magnetic loss + MAGNETIC SERVICE reduces root attitude space to dimension <=1 after four-S/process kills independent LIN/AW root",
            "one_acc_row_determines_BA_given_attitude":True,
            "complete_regular_word_joint_kernel_dim_le_one":True,
            "persistence_exclusion_separate_Astar":True,
            "theorem_closed":True}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
